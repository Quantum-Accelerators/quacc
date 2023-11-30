from ase.calculators.espresso import Espresso as _Espresso
from ase.calculators.espresso import EspressoProfile
from ase.calculators.espresso import EspressoTemplate as _EspressoTemplate

from quacc import SETTINGS
from quacc.calculators.espresso.io import read, write
from quacc.calculators.espresso.keys import ALL_KEYS
from quacc.calculators.espresso.utils import (construct_namelist,
                                              parse_pp_and_cutoff)
from quacc.utils.dicts import merge_dicts
from quacc.utils.files import load_yaml_calc


class EspressoTemplate(_EspressoTemplate):

    def __init__(self, binary):
        super().__init__()
        self.input_file = f'{binary}.in'
        self.output_file = f'{binary}.out'
        self.binary = binary

    def write_input(self, directory, atoms, parameters, properties):
        dst = directory / self.inputname
        write(dst,
              atoms,
              format=self.binary,
              properties=properties,
              **parameters)

    def read_results(self, directory):
            path = directory / self.outputname
            atoms = read(path, format=self.binary)
            return dict(atoms.calc.properties())


class Espresso(_Espresso):

    def __init__(self,
                 input_atoms = None,
                 preset = None,
                 template = None,
                 profile = None,
                 calc_defaults = None,
                 **kwargs):

        self.preset = preset
        self.input_atoms = input_atoms
        self.calc_defaults = calc_defaults
        #input_data = kwargs.pop('input_data', None)
        #pseudopotentials = kwargs.pop('pseudopotentials', None)
        #kpts = kwargs.pop('kpts', None)
        
        template = template or EspressoTemplate('pw')
        profile = profile or EspressoProfile(argv=
                                str(SETTINGS.ESPRESSO_CMD).split())
        
        kwargs = self._kwargs_handler(template.binary, **kwargs)
        
        super().__init__(
            profile = profile,
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