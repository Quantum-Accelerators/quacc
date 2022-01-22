from htase.util.calc import load_yaml_calc, set_magmoms
from htase.util.atoms import check_is_metal, get_highest_block
from ase.calculators.vasp import Vasp
from ase.calculators.vasp.setups import _setups_defaults as ase_default_setups
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.bandstructure import HighSymmKpath
import numpy as np
import os
import warnings
from pathlib import Path
from copy import deepcopy

DEFAULT_CALCS_DIR = os.path.join(
    Path(__file__).resolve().parent, "..", "..", "defaults", "user_calcs", "vasp"
)


def SmartVasp(
    atoms,
    preset=None,
    incar_copilot=True,
    force_gamma=True,
    copy_magmoms=True,
    mag_default=1.0,
    mag_cutoff=0.05,
    verbose=True,
    **kwargs,
):
    """
    This is a wrapper around the VASP calculator that adjusts INCAR parameters on-the-fly.
    Also supports several automatic k-point generation schemes from Pymatgen.

    Parameters
    ----------
    atoms : ase.Atoms
        The Atoms object to be used for the calculation.
    preset : str
        The path to a .yaml file containing a list of INCAR parameters to use as a "preset"
            for the calculator. If no filepath is present, it will look in htase/defaults/user_calcs, such
            that preset="bulk_base" is supported. It will append .yaml at the end if not present.
    incar_copilot : bool
        If True, the INCAR parameters will be adjusted if they go against the VASP manual.
            Defaults to True.
    force_gamma : bool
        If True, the automatic k-point generation schemes will default to gamma-centered.
            Defaults to True.
    copy_magmoms : bool
        If True, any pre-existing atoms.get_magnetic_moments() will be set in atoms.set_initial_magnetic_moments().
            Defaults to True. Set this to False if you want to use a preset's magnetic moments every time.
    mag_default : float
        Default magmom value for sites without one in the preset. Use 0.6 for MP settings or 1.0 for VASP default.
            Defaults to 1.0.
    mag_cutoff : None or float
        Set all initial magmoms to 0 if all have a magnitude below this value.
            Defaults to 0.05.
    verbose : bool
        If True, warnings will be raised when INCAR parameters are changed.
            Defaults to True.
    **kwargs :
        Additional arguments to be passed to the VASP calculator, e.g. xc='PBE', encut=520. Takes all valid
            ASE calculator arguments, in addition to those custom to HT-ASE.

    Returns
    -------
    atoms : ase.Atoms
        The ASE Atoms object with attached VASP calculator.
    """

    # Copy the atoms to a new object so we don't modify the original
    atoms = deepcopy(atoms)

    # Grab the pymatgen structure object in case we need it later
    struct = AseAtomsAdaptor().get_structure(atoms)

    # Get user-defined preset parameters for the calculator
    if preset:
        calc_preset = load_yaml_calc(preset, default_calcs_dir=DEFAULT_CALCS_DIR)[
            "inputs"
        ]
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
            user_calc_params["setups"], default_calcs_dir=DEFAULT_CALCS_DIR
        )["inputs"]["setups"]

    # If the user explicitly requests gamma = False, let's honor that
    # over force_gamma.
    if user_calc_params.get("gamma", None) is not None:
        user_gamma = user_calc_params["gamma"]
        if force_gamma is True:
            raise ValueError("force_gamma is True, but gamma is set to be False.")
    else:
        user_gamma = None

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
        kpts, gamma, reciprocal = _convert_auto_kpts(
            struct, auto_kpts, force_gamma=force_gamma
        )
        user_calc_params["kpts"] = kpts
        if reciprocal:
            user_calc_params["reciprocal"] = reciprocal
        if user_gamma is None:
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
    user_calc_params = _remove_unused_flags(user_calc_params)

    # Instantiate the calculator!
    calc = Vasp(**user_calc_params)

    # Handle INCAR swaps as needed
    if incar_copilot:
        calc = _calc_swaps(atoms, calc, auto_kpts=auto_kpts, verbose=verbose)

    # This is important! We want to make sure that setting
    # a new VASP parameter throws away the prior calculator results
    # otherwise we can't do things like run a new calculation
    # with atoms.get_potential_energy() after the transformation
    calc.discard_results_on_any_change = True

    # Set the calculator
    atoms.calc = calc

    return atoms


def _convert_auto_kpts(struct, auto_kpts, force_gamma=True):
    """
    Shortcuts for pymatgen k-point generation schemes.
    Options include: line_density (for band structures),
    reciprocal_density (by volume), grid_density (by number of atoms),
    max_mixed_density (max of reciprocal_density volume or atoms), and length_density (good for slabs).
    These are formatted as {"line_density": float}, {"reciprocal_density": float},
    {"grid_density": float}, {"vol_kkpa_density": [float, float]}, and
    {"length_density": [float, float, float]}.

    Args:
        struct (pymatgen.core.Structure): Pymatgen Structure object
        auto_kpts (str): String indicating the automatic k-point scheme
        force_gamma (bool): Force gamma-centered k-point grid
            Default: True

    Returns:
        kpts (List[Int]): ASE-compatible kpts argument
        gamma (bool): ASE gamma kwarg
        reciprocal (bool): ASE reciprocal kwarg

    """

    if auto_kpts.get("line_density", None):
        # TODO: Support methods other than latimer-munro
        kpath = HighSymmKpath(struct, path_type="latimer_munro")
        kpts, _ = kpath.get_kpoints(
            line_density=auto_kpts["line_density"], coords_are_cartesian=True
        )
        kpts = np.stack(kpts)
        reciprocal = True
        gamma = None

    else:
        reciprocal = None
        if auto_kpts.get("max_mixed_density", None):
            if len(auto_kpts["max_mixed_density"]) != 2:
                raise ValueError("Must specify two values for max_mixed_density.")

            if auto_kpts["max_mixed_density"][0] > auto_kpts["max_mixed_density"][1]:
                warnings.warn(
                    "Warning: It is not usual that kppvol > kppa. Please make sure you have chosen the right k-point densities.",
                    UserWarning,
                )
            pmg_kpts1 = Kpoints.automatic_density_by_vol(
                struct, auto_kpts["max_mixed_density"][0], force_gamma
            )
            pmg_kpts2 = Kpoints.automatic_density(
                struct, auto_kpts["max_mixed_density"][1], force_gamma
            )
            if np.product(pmg_kpts1.kpts[0]) >= np.product(pmg_kpts2.kpts[0]):
                pmg_kpts = pmg_kpts1
            else:
                pmg_kpts = pmg_kpts2
        elif auto_kpts.get("reciprocal_density", None):
            pmg_kpts = Kpoints.automatic_density_by_vol(
                struct, auto_kpts["reciprocal_density"], force_gamma
            )
        elif auto_kpts.get("grid_density", None):
            pmg_kpts = Kpoints.automatic_density(
                struct, auto_kpts["grid_density"], force_gamma
            )
        elif auto_kpts.get("length_density", None):
            if len(auto_kpts["length_density"]) != 3:
                raise ValueError("Must specify three values for length_density.")
            pmg_kpts = Kpoints.automatic_density_by_lengths(
                struct, auto_kpts["length_density"], force_gamma
            )
        else:
            raise ValueError(f"Unsupported k-point generation scheme: {auto_kpts}.")

        kpts = pmg_kpts.kpts[0]
        if pmg_kpts.style.name.lower() == "gamma":
            gamma = True
        else:
            gamma = False

    return kpts, gamma, reciprocal


def _remove_unused_flags(user_calc_params):
    """
    Removes unused flags in the INCAR, like EDIFFG if you are doing NSW = 0.

    Args:
        user_calc_params (dict): User-specified calculation parameters

    Returns:
        user_calc_params (dict): Adjusted user-specified calculation parameters
    """

    # Turn off opt flags if NSW = 0
    opt_flags = ("ediffg", "ibrion", "isif", "potim", "iopt")
    if user_calc_params.get("nsw", 0) == 0:
        for opt_flag in opt_flags:
            user_calc_params.pop(opt_flag, None)

    # Turn off +U flags if +U is not even used
    ldau_flags = (
        "ldau",
        "ldauu",
        "ldauj",
        "ldaul",
        "ldautype",
        "ldauprint",
        "ldau_luj",
    )
    if not user_calc_params.get("ldau", False) and not user_calc_params.get(
        "ldau_luj", None
    ):
        for ldau_flag in ldau_flags:
            user_calc_params.pop(ldau_flag, None)

    return user_calc_params


def _calc_swaps(atoms, calc, auto_kpts=None, verbose=True):
    """
    Swaps out bad INCAR flags.

    Args:
        atoms (ase.Atoms): Atoms object
        calc (ase.calculators.vasp.Vasp): ASE VASP calculator
        auto_kpts (bool): Whether to automatically set kpoints using
            one of the Pymatgen schemes.
            Default: None.
        verbose (bool): Whether to print out any time the input flags
            are adjusted.
            Default: True.

    Returns:
        calc (ase.calculators.vasp.Vasp): ASE VASP calculator with
            modified arguments.
    """
    is_metal = check_is_metal(atoms)
    max_block = get_highest_block(atoms)

    if (
        not calc.int_params["lmaxmix"] or calc.int_params["lmaxmix"] < 6
    ) and max_block == "f":
        if verbose:
            warnings.warn(
                "Copilot: Setting LMAXMIX = 6 because you have an f-element.",
                UserWarning,
            )
        calc.set(lmaxmix=6)
    elif (
        not calc.int_params["lmaxmix"] or calc.int_params["lmaxmix"] < 4
    ) and max_block == "d":
        if verbose:
            warnings.warn(
                "Copilot: Setting LMAXMIX = 4 because you have a d-element",
                UserWarning,
            )
        calc.set(lmaxmix=4)

    if (
        calc.bool_params["luse_vdw"]
        or calc.bool_params["lhfcalc"]
        or calc.bool_params["ldau"]
        or calc.dict_params["ldau_luj"]
        or calc.string_params["metagga"]
    ) and not calc.bool_params["lasph"]:
        if verbose:
            warnings.warn(
                "Copilot: Setting LASPH = True because you have a +U, vdW, meta-GGA, or hybrid calculation.",
                UserWarning,
            )
        calc.set(lasph=True)

    if (
        calc.bool_params["lasph"]
        and (not calc.int_params["lmaxtau"] or calc.int_params["lmaxtau"] < 8)
        and max_block == "f"
    ):
        if verbose:
            warnings.warn(
                "Copilot: Setting LMAXTAU = 8 because you have LASPH = True and an f-element.",
                UserWarning,
            )
        calc.set(lmaxtau=8)

    if calc.string_params["metagga"] and (
        not calc.string_params["algo"] or calc.string_params["algo"].lower() != "all"
    ):
        if verbose:
            warnings.warn(
                "Copilot: Setting ALGO = All because you have a meta-GGA calculation.",
                UserWarning,
            )
        calc.set(algo="all")

    if calc.bool_params["lhfcalc"] and (
        not calc.string_params["algo"]
        or calc.string_params["algo"].lower() not in ["all", "damped"]
    ):
        if is_metal:
            calc.set(algo="damped", time=0.5)
            if verbose:
                warnings.warn(
                    "Copilot: Setting ALGO = Damped, TIME = 0.5 because you have a hybrid calculation with a metal.",
                    UserWarning,
                )
        else:
            calc.set(algo="all")
            if verbose:
                warnings.warn(
                    "Copilot: Setting ALGO = All because you have a hybrid calculation.",
                    UserWarning,
                )

    if (
        is_metal
        and (calc.int_params["ismear"] and calc.int_params["ismear"] < 0)
        and (calc.int_params["nsw"] and calc.int_params["nsw"] > 0)
    ):
        if verbose:
            warnings.warn(
                "Copilot: You are relaxing a likely metal. Setting ISMEAR = 1 and SIGMA = 0.1.",
                UserWarning,
            )
        calc.set(ismear=1, sigma=0.1)

    if (
        calc.int_params["nedos"]
        and calc.int_params["ismear"] != -5
        and calc.int_params["nsw"] in (None, 0)
    ):
        if verbose:
            warnings.warn(
                "Copilot: Setting ISMEAR = -5 and SIGMA = 0.05 because you have a static DOS calculation.",
                UserWarning,
            )
        calc.set(ismear=-5, sigma=0.05)

    if calc.int_params["ismear"] == -5 and np.product(calc.kpts) < 4:
        if verbose:
            warnings.warn(
                "Copilot: Setting ISMEAR = 0 and SIGMA = 0.05 because you don't have enough k-points for ISMEAR = -5.",
                UserWarning,
            )
        calc.set(ismear=0, sigma=0.05)

    if (
        auto_kpts
        and auto_kpts.get("line_density", None)
        and (calc.int_params["ismear"] != 0 or calc.float_params["sigma"] > 0.01)
    ):
        if verbose:
            warnings.warn(
                "Copilot: Setting ISMEAR = 0 and SIGMA = 0.01 because you are doing a line mode calculation.",
                UserWarning,
            )
        calc.set(ismear=0, sigma=0.01)

    if (
        calc.float_params["kspacing"]
        and (calc.float_params["kspacing"] and calc.float_params["kspacing"] > 0.5)
        and calc.int_params["ismear"] == -5
    ):
        if verbose:
            warnings.warn(
                "Copilot: KSPACING is likely too large for ISMEAR = -5. Setting ISMEAR = 0 and SIGMA = 0.05.",
                UserWarning,
            )
        calc.set(ismear=0, sigma=0.05)
        pass

    if (
        calc.int_params["nsw"]
        and calc.int_params["nsw"] > 0
        and calc.bool_params["laechg"]
    ):
        if verbose:
            warnings.warn(
                "Copilot: Setting LAECHG = False because you have NSW > 0. LAECHG is not compatible with NSW > 0.",
                UserWarning,
            )
        calc.set(laechg=None)

    if calc.int_params["ldauprint"] in (None, 0) and (
        calc.bool_params["ldau"] or calc.dict_params["ldau_luj"]
    ):
        if verbose:
            warnings.warn(
                "Copilot: Setting LDAUPRINT = 1 because LDAU = True.", UserWarning
            )
        calc.set(ldauprint=1)

    if calc.special_params["lreal"] and calc.int_params["nsw"] in (None, 0, 1):
        if verbose:
            warnings.warn(
                "Copilot: Setting LREAL = False because you are running a static calculation. LREAL != False can be bad for energies.",
                UserWarning,
            )
        calc.set(lreal=False)

    if not calc.int_params["lorbit"] and (
        calc.int_params["ispin"] == 2
        or np.any(atoms.get_initial_magnetic_moments() != 0)
    ):
        if verbose:
            warnings.warn(
                "Copilot: Setting LORBIT = 11 because you have a spin-polarized calculation.",
                UserWarning,
            )
        calc.set(lorbit=11)

    if (
        (calc.int_params["ncore"] and calc.int_params["ncore"] > 1)
        or (calc.int_params["npar"] and calc.int_params["npar"] > 1)
    ) and (
        calc.bool_params["lhfcalc"] is True
        or calc.bool_params["lrpa"] is True
        or calc.bool_params["lepsilon"] is True
        or calc.int_params["ibrion"] in [5, 6, 7, 8]
    ):
        if verbose:
            warnings.warn(
                "Copilot: Setting NCORE = 1 because NCORE/NPAR is not compatible with this job type.",
                UserWarning,
            )
        calc.set(ncore=1)
        if calc.int_params.get("npar", None) is not None:
            calc.int_params["npar"] = None

    if (
        (calc.int_params["ncore"] and calc.int_params["ncore"] > 1)
        or (calc.int_params["npar"] and calc.int_params["npar"] > 1)
    ) and len(atoms) <= 4:
        if verbose:
            warnings.warn(
                "Copilot: Setting NCORE = 1 because you have a very small structure.",
                UserWarning,
            )
        calc.set(ncore=1)
        if calc.int_params.get("npar", None) is not None:
            calc.int_params["npar"] = None

    if calc.int_params["kpar"] and calc.int_params["kpar"] > np.prod(calc.kpts):
        if verbose:
            warnings.warn(
                "Copilot: Setting KPAR = 1 because you have too few k-points to parallelize.",
                UserWarning,
            )
        calc.set(kpar=1)

    if (
        calc.int_params["nsw"]
        and calc.int_params["nsw"] > 0
        and calc.int_params["isym"]
        and calc.int_params["isym"] > 0
    ):
        if verbose:
            warnings.warn(
                "Copilot: Setting ISYM = 0 because you are running a relaxation.",
                UserWarning,
            )
        calc.set(isym=0)

    if (
        calc.bool_params["lhfcalc"] is True
        and calc.int_params["isym"]
        and calc.int_params["isym"] in (1, 2)
    ):
        if verbose:
            warnings.warn(
                "Copilot: Setting ISYM = 3 because you are running a hybrid calculation.",
                UserWarning,
            )
        calc.set(isym=3)

    if calc.bool_params["luse_vdw"] and "ASE_VASP_VDW" not in os.environ:
        warnings.warn(
            "ASE_VASP_VDW was not set, yet you requested a vdW functional.",
            UserWarning,
        )

    return calc
