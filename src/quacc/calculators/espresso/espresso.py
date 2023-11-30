from ase.calculators.espresso import Espresso as _Espresso
from ase.calculators.espresso import EspressoTemplate as _EspressoTemplate
from ase.calculators.espresso import EspressoProfile, EspressoTemplate
from ase.io.espresso import construct_namelist

from quacc import SETTINGS
from quacc.calculators.espresso.utils import parse_pp_and_cutoff
from quacc.utils.dicts import merge_dicts
from quacc.utils.files import load_yaml_calc
from quacc.calculators.espresso.io import write, read


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
                 template = _EspressoTemplate(),
                 **kwargs):

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
    
        # By default we fall back on ase.espresso.EspressoTemplate
        self.template = template

    def _kwargs_handler(self, **kwargs):
        kwargs['input_data'] = construct_namelist(
            kwargs.get('input_data', None))
        if self.preset:
            config = load_yaml_calc(
                SETTINGS.ESPRESSO_PRESET_DIR / f"{self.preset}"
            )
            preset_pp = parse_pp_and_cutoff(config, self.input_atoms)
            return merge_dicts(preset_pp, kwargs)