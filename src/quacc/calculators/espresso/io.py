"""I/O utilities for the espresso calculator."""

# Three ways to call a function based on a string:

# 1. use the globals() dictionary
# 2. use sys.modules[__name__] and getattr
# 3. keep track of all functions in a dictionary

def write(filename,
          atoms,
          format = 'espresso-in',
          properties = None,
          **kwargs):
    pass
    #lines = write_namelist(filename, parameters)
    #specific_lines = func_dict[format][1](**kwargs)
    #lines += specific_lines
    #with open(filename, 'w') as fd:
    #    fd.write(''.join(lines))

def read(filename,
         format = 'espresso-out'):
    pass
    #with open(filename, 'r') as fd:
    #    lines = fd.readlines()
    #return func_dict[format][0](lines, **kwargs)

def write_namelist(parameters):
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
                pwi.append('   {0:16} = {1!r:}\n'.format(key, value))
        pwi.append('/\n')  # terminate section
    pwi.append('\n')
    return pwi

def write_ph_specifics():
    return ''

def read_ph_specifics():
    return ''

func_dict = {
    'ph': (read_ph_specifics, write_ph_specifics),
}