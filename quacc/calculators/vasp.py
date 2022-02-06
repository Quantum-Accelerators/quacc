import os
from copy import deepcopy
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from ase.atoms import Atoms
from ase.calculators.vasp import Vasp
from ase.calculators.vasp.setups import _setups_defaults as ase_default_setups
from ase.constraints import FixAtoms
from quacc.util.yaml import load_yaml_calc
from quacc.util.atoms import set_magmoms
from quacc.defaults.calcs import vasp as vasp_defaults
from quacc.calculators.vasp_utils import (
    manage_environment,
    convert_auto_kpts,
    remove_unused_flags,
    calc_swaps,
)

DEFAULT_CALCS_DIR = os.path.dirname(vasp_defaults.__file__)


def SmartVasp(
    atoms: Atoms,
    custodian: bool = True,
    preset: None | str = None,
    incar_copilot: bool = True,
    copy_magmoms: bool = True,
    mag_default: float = 1.0,
    mag_cutoff: None | float = 0.05,
    verbose: bool = True,
    **kwargs,
) -> Atoms:
    """
    This is a wrapper around the VASP calculator that adjusts INCAR parameters on-the-fly.
    Also supports several automatic k-point generation schemes from Pymatgen.

    Parameters
    ----------
    atoms
        The Atoms object to be used for the calculation.
    custodian
        Whether to use Custodian to run VASP. If True, the Custodian settings will (by default) be
        adopted from htase.defaults.custodian_settings.vasp_custodian_settings.yaml. To override
        these default Custodian settings, you can define your own .yaml file and set the path in
        the VASP_CUSTODIAN_SETTINGS environment variable at runtime. If set to False, ASE will
        run VASP without Custodian, in which case you must have a valid ASE_VASP_COMMAND set per
        the ASE documentation.
    preset
        The path to a .yaml file containing a list of INCAR parameters to use as a "preset"
        for the calculator. If no filepath is present, it will look in htase/defaults/user_calcs, such
        that preset="BulkRelaxSet" is supported. It will append .yaml at the end if not present.
        Note that any specific kwargs take precedence over the flags set in the preset dictionary.
    incar_copilot
        If True, the INCAR parameters will be adjusted if they go against the VASP manual.
    copy_magmoms
        If True, any pre-existing atoms.get_magnetic_moments() will be set in atoms.set_initial_magnetic_moments().
        Set this to False if you want to use a preset's magnetic moments every time.
    mag_default
        Default magmom value for sites without one in the preset. Use 0.6 for MP settings or 1.0 for VASP default.
    mag_cutoff
        Set all initial magmoms to 0 if all have a magnitude below this value.
    verbose
        If True, warnings will be raised when INCAR parameters are changed.
    **kwargs
        Additional arguments to be passed to the VASP calculator, e.g. xc='PBE', encut=520. Takes all valid
        ASE calculator arguments, in addition to those custom to HT-ASE.

    Returns
    -------
    Atoms
        The ASE Atoms object with attached VASP calculator.
    """

    # Copy the atoms to a new object so we don't modify the original
    atoms = deepcopy(atoms)

    # Grab the pymatgen structure object in case we need it later
    struct = AseAtomsAdaptor().get_structure(atoms)

    # Check constraints
    if custodian and atoms.constraints:
        if ~np.all([isinstance(c, FixAtoms) for c in atoms.constraints]):
            raise ValueError(
                "Atoms object has a constraint that is not compatible with Custodian"
            )

    # Get VASP executable command, if necessary, and specify child environment
    # variables
    command = manage_environment(custodian)

    # Get user-defined preset parameters for the calculator
    if preset:
        _, ext = os.path.splitext(preset)
        if not ext:
            preset += ".yaml"
        calc_preset = load_yaml_calc(os.path.join(DEFAULT_CALCS_DIR, preset))["inputs"]
    else:
        calc_preset = {}

    # Collect all the calculator parameters and prioritize the kwargs
    # in the case of duplicates.
    user_calc_params = {**calc_preset, **kwargs}
    none_keys = [k for k, v in user_calc_params.items() if v is None]
    for none_key in none_keys:
        del user_calc_params[none_key]

    # Allow the user to use setups='mysetups.yaml' to load in a custom setups
    # from a YAML file
    if (
        isinstance(user_calc_params.get("setups", None), str)
        and user_calc_params["setups"] not in ase_default_setups
    ):
        user_calc_params["setups"] = load_yaml_calc(
            os.path.join(DEFAULT_CALCS_DIR, user_calc_params["setups"])
        )["inputs"]["setups"]

    # If the preset has auto_kpts but the user explicitly requests kpts, then
    # we should honor that.
    if kwargs.get("kpts", None) and calc_preset.get("auto_kpts", None):
        del user_calc_params["auto_kpts"]

    # If the preset has ediff_per_atom but the user explicitly requests ediff, then
    # we should honor that.
    if kwargs.get("ediff", None) and calc_preset.get("ediff_per_atom", None):
        del user_calc_params["ediff_per_atom"]

    # Handle special arguments in the user calc parameters that
    # ASE does not natively support
    if user_calc_params.get("elemental_magmoms", None) is not None:
        elemental_mags_dict = user_calc_params["elemental_magmoms"]
        del user_calc_params["elemental_magmoms"]
    else:
        elemental_mags_dict = None
    if user_calc_params.get("auto_kpts", None) is not None:
        auto_kpts = user_calc_params["auto_kpts"]
        del user_calc_params["auto_kpts"]
    else:
        auto_kpts = None
    if user_calc_params.get("auto_dipole", None) is not None:
        auto_dipole = user_calc_params["auto_dipole"]
        del user_calc_params["auto_dipole"]
    else:
        auto_dipole = None
    if user_calc_params.get("ediff_per_atom", None) is not None:
        ediff_per_atom = user_calc_params["ediff_per_atom"]
        del user_calc_params["ediff_per_atom"]
    else:
        ediff_per_atom = None

    # Make automatic k-point mesh
    if auto_kpts:
        kpts, gamma, reciprocal = convert_auto_kpts(struct, auto_kpts)
        user_calc_params["kpts"] = kpts
        if reciprocal and user_calc_params.get("reciprocal", None) is None:
            user_calc_params["reciprocal"] = reciprocal
        if user_calc_params.get("gamma", None) is None:
            user_calc_params["gamma"] = gamma

    # Add dipole corrections if requested
    if auto_dipole:
        com = atoms.get_center_of_mass(scaled=True)
        if "dipol" not in user_calc_params:
            user_calc_params["dipol"] = com
        if "idipol" not in user_calc_params:
            user_calc_params["idipol"] = 3
        if "ldipol" not in user_calc_params:
            user_calc_params["ldipol"] = True

    # Handle ediff_per_atom if present
    if ediff_per_atom:
        user_calc_params["ediff"] = ediff_per_atom * len(atoms)

    # Set magnetic moments
    atoms = set_magmoms(
        atoms,
        elemental_mags_dict=elemental_mags_dict,
        copy_magmoms=copy_magmoms,
        mag_default=mag_default,
        mag_cutoff=mag_cutoff,
    )
    if ~np.all(
        [isinstance(m, (int, float)) for m in atoms.get_initial_magnetic_moments()]
    ):
        raise ValueError(
            "Magnetic moments must be specified as a list of floats or ints.",
            atoms.get_initial_magnetic_moments().tolist(),
        )

    # Remove unused INCAR flags
    user_calc_params = remove_unused_flags(user_calc_params)

    # Instantiate the calculator!
    calc = Vasp(command=command, **user_calc_params)

    # Handle INCAR swaps as needed
    if incar_copilot:
        calc = calc_swaps(atoms, calc, auto_kpts=auto_kpts, verbose=verbose)

    # This is important! We want to make sure that setting
    # a new VASP parameter throws away the prior calculator results
    # otherwise we can't do things like run a new calculation
    # with atoms.get_potential_energy() after the transformation
    calc.discard_results_on_any_change = True

    # Set the calculator
    atoms.calc = calc

    return atoms
