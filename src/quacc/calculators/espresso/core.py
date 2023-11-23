from quacc.calculators.espresso.utils import construct_namelist
from quacc.calculators.espresso.io import _write, _read
from quacc.utils.files import load_yaml_calc

# Very simple class that can run calculations ASE style
# I don't think we even need atoms object here.
class EspressoCalculator:
    # These are not supposed to change
    def __init__(self):
        self.input_file = f'{self._BIN}.in'
        self.output_file = f'{self._BIN}.out'
        self._calc_defaults = self.to_namelist(self._calc_defaults)

    def write_input(self, parameters, **kwargs):
        _write(self.input_file, parameters, format = self._BIN, **kwargs)

    def run(self):
        from os import environ
        from subprocess import check_call
        argv = list(self._CMD) + ['-in', self.input_file]
        with open(self.output_file, 'wb') as fd:
            check_call(argv, stdout=fd, env=environ)
    
    def execute(self, preset = None, **kwargs):
        parameters = self.input_handler(self.input_data, preset)
        self.write_input(parameters, **kwargs)
        self.run()

    def read_results(self, **kwargs):
        return _read(self.output_file, format = self._BIN, **kwargs)
    
    # Could be an abstract method
    def input_handler(self, parameters, preset):
        parameters = construct_namelist(parameters, keys=self._KEYS)
        parameters = self._calc_defaults | parameters
        if preset:
            parameters = self.load_presets(preset) | parameters
        return parameters
    
    def load_presets(self, preset):
        parameters = load_yaml_calc(preset)
        return parameters
    
    def to_namelist(self, parameters, keys = None):
        keys = self._KEYS if keys is None else keys
        return construct_namelist(parameters, keys=keys)


