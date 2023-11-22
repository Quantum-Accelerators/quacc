"""I/O utilities for the espresso calculator."""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Any


# Three ways to call a function based on a string:

# 1. use the globals() dictionary
# 2. use sys.modules[__name__] and getattr
# 3. keep track of all functions in a dictionary

def _write(filename, parameters, *, format, **kwargs):
    lines = write_namelist(filename, parameters)

    with open(filename, 'w') as fd:
        fd.write(''.join(lines))

def write_namelist(filename: str | Path, parameters: dict):
    # This function aim to write basic input keywords to a file
    # For additional cards (binary specifics) this is down in the child class
    # Stolen from ase.io.espresso, maybe this could be made to a function "write_namelist"
    # and imported here? See with the ASE gods.
    pwi = []
    for section in parameters:
        pwi.append('&{0}\n'.format(section.upper()))
        for key, value in parameters[section].items():
            if value is True:
                pwi.append('   {0:16} = .true.\n'.format(key))
            elif value is False:
                pwi.append('   {0:16} = .false.\n'.format(key))
            else:
                # repr format to get quotes around strings
                pwi.append('   {0:16} = {1!r:}\n'.format(key, value))
        pwi.append('/\n')  # terminate section
    pwi.append('\n')
    return pwi

def write_ph_specifics():
    pass