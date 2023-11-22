from quacc import SETTINGS
from quacc.calculators.espresso.core import EspressoCalculator
from quacc.calculators.espresso.keys import (DYNMAT_KEYS, NEB_KEYS, PH_KEYS,
                                             Q2R_KEYS, MATDYN_KEYS)

# Mini classes that might be useful for binaries
# that have very specific requirements

class Ph(EspressoCalculator):
    # These are not supposed to change
    _BIN = 'ph'
    _KEYS = PH_KEYS
    _CMD = SETTINGS['ESPRESSO_PH_CMD']

    def __init__(self, input_data):
        self.input_data = input_data
        super().__init__()

class Neb(EspressoCalculator):

    _BIN = 'neb'
    _KEYS = NEB_KEYS
    _CMD = SETTINGS['ESPRESSO_NEB_CMD']

    def __init__(self, input_data):
        self.input_data = input_data
        super().__init__()
        self.input_file = 'neb.dat'

class Q2r(EspressoCalculator):
    # These are not supposed to change
    _BIN = 'q2r'
    _KEYS = Q2R_KEYS
    _CMD = SETTINGS['ESPRESSO_Q2R_CMD']

    def __init__(self, input_data):
        self.input_data = input_data
        super().__init__()

class Dynmat(EspressoCalculator):
    # These are not supposed to change
    _BIN = 'dynmat'
    _KEYS = DYNMAT_KEYS
    _CMD = SETTINGS['ESPRESSO_DYNMAT_CMD']

    def __init__(self, input_data):
        self.input_data = input_data
        super().__init__()

class Matdyn(EspressoCalculator):
    # These are not supposed to change
    _BIN = 'matdyn'
    _KEYS = MATDYN_KEYS
    _CMD = SETTINGS['ESPRESSO_MATDYN_CMD']

    def __init__(self, input_data):
        self.input_data = input_data
        super().__init__()