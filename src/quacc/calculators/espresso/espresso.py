from ase.calculators.espresso import Espresso as _Espresso
from ase.calculators.espresso import EspressoProfile
from ase.io.espresso import construct_namelist

from quacc import SETTINGS
from quacc.calculators.espresso.utils import parse_pp_and_cutoff
from quacc.utils.dicts import merge_dicts
from quacc.utils.files import load_yaml_calc


class Espresso(_Espresso):

    def __init__(self, input_atoms = None, preset = None, **kwargs):

        self.preset = preset
        self.input_atoms = input_atoms

        kwargs = self._kwargs_handler(**kwargs)

        input_data = kwargs.pop('input_data', None)
        profile = kwargs.pop('profile',
                             EspressoProfile(argv=
                             str(SETTINGS.ESPRESSO_CMD).split()))
        pseudopotentials = kwargs.pop('pseudopotentials', None)
        kpts = kwargs.pop('kpts', None)

        super().__init__(
            profile = profile,
            input_data = input_data,
            pseudopotentials = pseudopotentials,
            kpts = kpts, 
            **kwargs)

    def _kwargs_handler(self, **kwargs):
        kwargs['input_data'] = construct_namelist(
            kwargs.get('input_data', None))
        if self.preset:
            config = load_yaml_calc(
                SETTINGS.ESPRESSO_PRESET_DIR / f"{self.preset}"
            )
            preset_pp = parse_pp_and_cutoff(config, self.input_atoms)
            return merge_dicts(preset_pp, kwargs)