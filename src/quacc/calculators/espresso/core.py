from quacc.calculators.espresso.utils import construct_namelist
from quacc.calculators.espresso.io import _write, _read

# Very simple class that can run calculations ASE style
# I don't think we even need atoms object here.
class EspressoCalculator:
    # These are not supposed to change
    def __init__(self):
        self.input_file = f'{self._BIN}.in'
        self.output_file = f'{self._BIN}.out'

    def write_input(self, parameters, **kwargs):
        _write(self.input_file, parameters, format = self._BIN, **kwargs)

    def run(self):
        from os import environ
        from subprocess import check_call
        argv = list(self._CMD) + ['-in', self.input_file]
        with open(self.output_file, 'wb') as fd:
            check_call(argv, stdout=fd, env=environ)
    
    def execute(self, **kwargs):
        parameters = self.to_namelist(self.input_data)
        self.write_input(parameters, **kwargs)
        self.run()

    def read_results(self, **kwargs):
        return _read(self.output_file, format = self._BIN, **kwargs)

    def to_namelist(self):
        return construct_namelist(self.input_data, keys=self._KEYS)