from typing import TYPE_CHECKING

from ase import Atoms
from ase.calculators.espresso import Espresso as Espresso_
from ase.calculators.espresso import EspressoProfile
from ase.calculators.espresso import EspressoTemplate as EspressoTemplate_

from quacc import SETTINGS
from quacc.calculators.espresso.io import read, write
from quacc.calculators.espresso.keys import ALL_KEYS
from quacc.calculators.espresso.utils import (construct_namelist,
                                              parse_pp_and_cutoff)
from quacc.utils.dicts import merge_dicts
from quacc.utils.files import load_yaml_calc

if TYPE_CHECKING:
    from typing import Any


class EspressoTemplate(EspressoTemplate_):

    def __init__(self,
                 binary: str = 'pw'):
        super().__init__()
        self.inputname = f'{binary}.in'
        self.outputname = f'{binary}.out'
        self.binary = binary

    def write_input(self, profile, directory,
                    atoms, parameters, properties):
        dst = directory / self.inputname
        write(dst,
              atoms,
              format=self.binary,
              properties=properties,
              pseudo_dir=str(profile.pseudo_path),
              **parameters)

    def read_results(self, directory):
        path = directory / self.outputname
        atoms = read(path, format=self.binary)
        return dict(atoms.calc.properties())


class Espresso(Espresso_):

    def __init__(self,
                 input_atoms: Atoms = None,
                 preset: str = None,
                 template: EspressoTemplate = None,
                 profile: EspressoProfile = None,
                 calc_defaults: dict[str | Any] = None,
                 parallel_info: dict[str | Any] = None,
                 **kwargs):

        self.preset = preset
        self.input_atoms = input_atoms
        self.calc_defaults = calc_defaults

        template = template or EspressoTemplate('pw')

        kwargs = self._kwargs_handler(template.binary, **kwargs)

        pseudo_path = kwargs['input_data']['control'].get(
            'pseudo_dir', str(SETTINGS.ESPRESSO_PP_PATH))

        profile = profile or EspressoProfile(
            binary=str(SETTINGS.ESPRESSO_CMD),
            parallel_info=parallel_info,
            pseudo_path=pseudo_path)

        super().__init__(
            profile=profile,
            parallel_info=parallel_info,
            **kwargs)

        self.template = template

    def _kwargs_handler(self, binary, **kwargs):
        keys = ALL_KEYS[binary]
        kwargs['input_data'] = construct_namelist(
            kwargs.get('input_data', None), keys=keys)
        self.calc_defaults['input_data'] = construct_namelist(
            self.calc_defaults['input_data'], keys=keys)
        # Would be nice to change the merge_dict function so that
        # it is fully compatible with the Namelist class. I believe
        # changing 'dict or {}' would do.
        if self.preset:
            config = load_yaml_calc(
                SETTINGS.ESPRESSO_PRESET_DIR / f"{self.preset}"
            )
            preset_pp = parse_pp_and_cutoff(config, self.input_atoms)
            kwargs = merge_dicts(preset_pp, kwargs)
        kwargs = merge_dicts(self.calc_defaults, kwargs)
        return kwargs
