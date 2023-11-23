from quacc import SETTINGS
from quacc.calculators.espresso.core import EspressoCalculator
from quacc.calculators.espresso.keys import (DYNMAT_KEYS, NEB_KEYS, PH_KEYS,
                                             Q2R_KEYS, MATDYN_KEYS, PW_KEYS)
from quacc.calculators.espresso.utils import parse_pw_preset

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

    # I think it makes sense to have the defaults
    # here. They are not supposed to change.
    # Users can have a look at them by doing:
    # from quacc.calculators.espresso.calculators import Neb
    # Neb._calc_defaults without looking at the source code

    # Note to self: 'When a child class method is created with the same name and signature as one in the parent, the child's method takes precedence' didn't know that, praise python
    _BIN = 'neb'
    _KEYS = NEB_KEYS
    _CMD = SETTINGS['ESPRESSO_NEB_CMD']

    _calc_defaults = {
        'input_data': {
            'path': {
                'string_method': 'neb',
                'restart_mode': 'from_scratch',
                'opt_scheme': 'broyden',
                'CI_scheme': 'no-CI',
                'first_last_opt': False,
                "num_of_images": 6
            }
        }
    }

    def __init__(self, neb_input_data, pw_input_data):
        self.input_data = neb_input_data
        self.additional_data = pw_input_data
        super().__init__()
        self.input_file = 'neb.dat'

    # It would be nice to be able to run sanity_checks BEFORE
    # the calculation is submitted to any kind of HPC
    # This could also raise warnings to remind the user
    # about known problems i.e. if first_last_opt is not True
    # and the first and last images contains spin state. Running
    # the first singlepoint with fresh wavefunctions might end up
    # in a different spin state, screwing up the whole NEB calculation.

    # Probably the best way is to have a recipe that only does sanity checks
    # with a different runner (local) and then the actual calculation with
    # the HPC runner. This way we can also check if the calculation is
    # going to fail.

    # Ideally quacc should completely prevent the first kind of error
    #   1 - "input file" error
    # quacc will should be able to partially help the second kind of error
    #   2 - "input file" is correct but parallelization or hpc configuration
    #      is wrong i.e. if image parallelization is used but the number of
    #      images is not a multiple of the number of cores or similar
    # quacc will not be able to prevent the third kind of error
    #   3 - everything is correct but the calculation fail for non-technical reasons
    # Although quacc should be able to detect it and fail the calculation. It is also
    # necessary to be able to restart the calculation in the simplest and best way authorized by each binary.
    # i.e. if the calculation fails because of a walltime limit, quacc should be able to
    # restart the calculation without too much interaction from the user.
    def sanity_checks(self):
        if len(self.atoms_list) < 2:
            raise ValueError(
                'You need at least 2 images to run a NEB calculation.')
    
    def input_handler(self, parameters, preset):
        parameters = self.to_namelist(parameters)
        parameters = self._calc_defaults | parameters
        if preset:
            preset_conf = self.load_presets(preset)
            preset_conf = parse_pw_preset(preset_conf)
        else:
            preset_conf = {}
        for i in self.additional_data.keys():
            input_data = self.to_namelist(
                self.additional_data[i]['input_data'],
                keys=PW_KEYS)
            self.additional_data[i]['input_data'] = input_data
            self.additional_data[i] = preset_conf | self.additional_data[i]

        

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