import inspect
import os
import warnings
from typing import Any, Dict, List, Tuple

import numpy as np
from ase.atoms import Atoms
from ase.calculators.vasp import Vasp
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.symmetry.bandstructure import HighSymmKpath

from quacc.custodian import vasp as custodian_vasp
from quacc.defaults import custodian_settings
from quacc.util.atoms import check_is_metal, get_highest_block


def manage_environment(custodian: bool = True) -> str:
    """
    Manage the environment for the VASP calculator.

    Parameters
    ----------
    custodian
        If True, Custodian will be used to run VASP.

    Returns
    -------
    str
        The command flag to pass to the Vasp calculator.
    """

    # Check ASE environment variables
    if "VASP_PP_PATH" not in os.environ:
        warnings.warn(
            "The VASP_PP_PATH environment variable must point to the library of VASP pseudopotentials. See the ASE Vasp calculator documentation for details.",
        )

    # Check if Custodian should be used and confirm environment variables are set
    if custodian:
        if "VASP_CUSTODIAN_SETTINGS" in os.environ:
            custodian_yaml = os.environ["VASP_CUSTODIAN_SETTINGS"]
        else:
            custodian_yaml = os.path.join(
                os.path.dirname(custodian_settings.__file__),
                "vasp_custodian_settings.yaml",
            )
            os.environ["VASP_CUSTODIAN_SETTINGS"] = custodian_yaml
        if not os.path.isfile(custodian_yaml):
            raise FileNotFoundError(f"{custodian_yaml} not found.")

        # Return the command flag
        run_vasp_custodian_file = os.path.abspath(inspect.getfile(custodian_vasp))
        command = f"python {run_vasp_custodian_file}"
    else:
        if "ASE_VASP_COMMAND" not in os.environ and "VASP_SCRIPT" not in os.environ:
            warnings.warn(
                "ASE_VASP_COMMAND or VASP_SCRIPT must be set in the environment to run VASP. See the ASE Vasp calculator documentation for details."
            )
        command = None

    return command


def convert_auto_kpts(
    atoms: Atoms,
    auto_kpts: None
    | Dict[str, float]
    | Dict[str, List[Tuple[float, float]]]
    | Dict[str, List[Tuple[float, float, float]]],
    force_gamma: bool = True,
) -> Tuple[List[Tuple[int, int, int]], None | bool, None | bool]:
    """
    Shortcuts for pymatgen k-point generation schemes.
    Options include: line_density (for band structures),
    reciprocal_density (by volume), grid_density (by number of atoms),
    max_mixed_density (max of reciprocal_density volume or atoms), and length_density (good for slabs).
    These are formatted as {"line_density": float}, {"reciprocal_density": float},
    {"grid_density": float}, {"vol_kkpa_density": [float, float]}, and
    {"length_density": [float, float, float]}.

    Parameters
    ----------
    .Atoms
        ASE Atoms object
    auto_kpts
        Dictionary describing the automatic k-point scheme
    force_gamma
        Whether a gamma-centered mesh should be returned

    Returns
    -------
    List[int, int, int]
        List of k-points for use with the ASE Vasp calculator
    Optional[bool]
        The gamma command for use with the ASE Vasp calculator
    Optional[bool]
        The reciprocal command for use with the ASE Vasp calculator
    """
    struct = AseAtomsAdaptor.get_structure(atoms)

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
                )
            pmg_kpts1 = Kpoints.automatic_density_by_vol(
                struct, auto_kpts["max_mixed_density"][0], force_gamma=force_gamma
            )
            pmg_kpts2 = Kpoints.automatic_density(
                struct, auto_kpts["max_mixed_density"][1], force_gamma=force_gamma
            )
            if np.product(pmg_kpts1.kpts[0]) >= np.product(pmg_kpts2.kpts[0]):
                pmg_kpts = pmg_kpts1
            else:
                pmg_kpts = pmg_kpts2
        elif auto_kpts.get("reciprocal_density", None):
            pmg_kpts = Kpoints.automatic_density_by_vol(
                struct, auto_kpts["reciprocal_density"], force_gamma=force_gamma
            )
        elif auto_kpts.get("grid_density", None):
            pmg_kpts = Kpoints.automatic_density(
                struct, auto_kpts["grid_density"], force_gamma=force_gamma
            )
        elif auto_kpts.get("length_density", None):
            if len(auto_kpts["length_density"]) != 3:
                raise ValueError("Must specify three values for length_density.")
            pmg_kpts = Kpoints.automatic_density_by_lengths(
                struct, auto_kpts["length_density"], force_gamma=force_gamma
            )
        else:
            raise ValueError(f"Unsupported k-point generation scheme: {auto_kpts}.")

        kpts = pmg_kpts.kpts[0]
        gamma = pmg_kpts.style.name.lower() == "gamma"

    return kpts, gamma, reciprocal


def remove_unused_flags(user_calc_params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Removes unused flags in the INCAR, like EDIFFG if you are doing NSW = 0.

    Parameters
    ----------
    user_calc_params
        User-specified calculation parameters

    Returns
    -------
    Dict
        Adjusted user-specified calculation parameters
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


def calc_swaps(
    atoms: Atoms,
    calc: Vasp,
    auto_kpts: None
    | Dict[str, float]
    | Dict[str, List[Tuple[float, float]]]
    | Dict[str, List[Tuple[float, float, float]]],
    verbose: bool = True,
) -> Vasp:
    """
    Swaps out bad INCAR flags.

    Parameters
    ----------
    .Atoms
        Atoms object
    .Vasp
        ASE VASP calculator
    auto_kpts
        The automatic k-point scheme dictionary
    verbose
        Whether to print out any time the input flags are adjusted.

    Returns
    -------
    .Vasp
        ASE VASP calculator with modified arguments.
    """
    is_metal = check_is_metal(atoms)
    max_block = get_highest_block(atoms)

    if (
        not calc.int_params["lmaxmix"] or calc.int_params["lmaxmix"] < 6
    ) and max_block == "f":
        if verbose:
            warnings.warn("Copilot: Setting LMAXMIX = 6 because you have an f-element.")
        calc.set(lmaxmix=6)
    elif (
        not calc.int_params["lmaxmix"] or calc.int_params["lmaxmix"] < 4
    ) and max_block == "d":
        if verbose:
            warnings.warn("Copilot: Setting LMAXMIX = 4 because you have a d-element")
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
                "Copilot: Setting LASPH = True because you have a +U, vdW, meta-GGA, or hybrid calculation."
            )
        calc.set(lasph=True)

    if (
        calc.bool_params["lasph"]
        and (not calc.int_params["lmaxtau"] or calc.int_params["lmaxtau"] < 8)
        and max_block == "f"
    ):
        if verbose:
            warnings.warn(
                "Copilot: Setting LMAXTAU = 8 because you have LASPH = True and an f-element."
            )
        calc.set(lmaxtau=8)

    if calc.string_params["metagga"] and (
        not calc.string_params["algo"] or calc.string_params["algo"].lower() != "all"
    ):
        if verbose:
            warnings.warn(
                "Copilot: Setting ALGO = All because you have a meta-GGA calculation."
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
                    "Copilot: Setting ALGO = Damped, TIME = 0.5 because you have a hybrid calculation with a metal."
                )
        else:
            calc.set(algo="all")
            if verbose:
                warnings.warn(
                    "Copilot: Setting ALGO = All because you have a hybrid calculation."
                )

    if (
        is_metal
        and (calc.int_params["ismear"] and calc.int_params["ismear"] < 0)
        and (calc.int_params["nsw"] and calc.int_params["nsw"] > 0)
    ):
        if verbose:
            warnings.warn(
                "Copilot: You are relaxing a likely metal. Setting ISMEAR = 1 and SIGMA = 0.1."
            )
        calc.set(ismear=1, sigma=0.1)

    if (
        calc.int_params["nedos"]
        and calc.int_params["ismear"] != -5
        and calc.int_params["nsw"] in (None, 0)
    ):
        if verbose:
            warnings.warn(
                "Copilot: Setting ISMEAR = -5 because you have a static DOS calculation."
            )
        calc.set(ismear=-5)
        if calc.float_params.get("sigma", 0.2) > 0.05:
            if verbose:
                    warnings.warn(
                        "Copilot: Setting SIGMA = 0.05 because a static DOS was requested with SIGMA > 0.05."
                    )
            calc.set(sigma=0.05)
            

    if (
        calc.int_params["ismear"] == -5
        and np.product(calc.kpts) < 4
        and calc.float_params["kspacing"] is None
    ):
        if verbose:
            warnings.warn(
                "Copilot: Setting ISMEAR = 0 and SIGMA = 0.05 because you don't have enough k-points for ISMEAR = -5."
            )
        calc.set(ismear=0, sigma=0.05)

    if (
        auto_kpts
        and auto_kpts.get("line_density", None)
        and (calc.int_params["ismear"] != 0 or calc.float_params["sigma"] > 0.01)
    ):
        if verbose:
            warnings.warn(
                "Copilot: Setting ISMEAR = 0 and SIGMA = 0.01 because you are doing a line mode calculation."
            )
        calc.set(ismear=0, sigma=0.01)

    if (
        calc.float_params["kspacing"]
        and (calc.float_params["kspacing"] and calc.float_params["kspacing"] > 0.5)
        and calc.int_params["ismear"] == -5
    ):
        if verbose:
            warnings.warn(
                "Copilot: KSPACING is likely too large for ISMEAR = -5. Setting ISMEAR = 0 and SIGMA = 0.05."
            )
        calc.set(ismear=0, sigma=0.05)

    if (
        calc.int_params["nsw"]
        and calc.int_params["nsw"] > 0
        and calc.bool_params["laechg"]
    ):
        if verbose:
            warnings.warn(
                "Copilot: Setting LAECHG = False because you have NSW > 0. LAECHG is not compatible with NSW > 0."
            )
        calc.set(laechg=None)

    if calc.int_params["ldauprint"] in (None, 0) and (
        calc.bool_params["ldau"] or calc.dict_params["ldau_luj"]
    ):
        if verbose:
            warnings.warn("Copilot: Setting LDAUPRINT = 1 because LDAU = True.")
        calc.set(ldauprint=1)

    if calc.special_params["lreal"] and calc.int_params["nsw"] in (None, 0, 1):
        if verbose:
            warnings.warn(
                "Copilot: Setting LREAL = False because you are running a static calculation. LREAL != False can be bad for energies."
            )
        calc.set(lreal=False)

    if not calc.int_params["lorbit"] and (
        calc.int_params["ispin"] == 2
        or np.any(atoms.get_initial_magnetic_moments() != 0)
    ):
        if verbose:
            warnings.warn(
                "Copilot: Setting LORBIT = 11 because you have a spin-polarized calculation."
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
                "Copilot: Setting NCORE = 1 because NCORE/NPAR is not compatible with this job type."
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
                "Copilot: Setting NCORE = 1 because you have a very small structure."
            )
        calc.set(ncore=1)
        if calc.int_params.get("npar", None) is not None:
            calc.int_params["npar"] = None

    if (
        calc.int_params["kpar"]
        and calc.int_params["kpar"] > np.prod(calc.kpts)
        and calc.float_params["kspacing"] is None
    ):
        if verbose:
            warnings.warn(
                "Copilot: Setting KPAR = 1 because you have too few k-points to parallelize."
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
                "Copilot: Setting ISYM = 0 because you are running a relaxation."
            )
        calc.set(isym=0)

    if calc.bool_params["lhfcalc"] is True and calc.int_params["isym"] in (1, 2):
        if verbose:
            warnings.warn(
                "Copilot: Setting ISYM = 3 because you are running a hybrid calculation."
            )
        calc.set(isym=3)

    if calc.bool_params["luse_vdw"] and "ASE_VASP_VDW" not in os.environ:
        warnings.warn("ASE_VASP_VDW was not set, yet you requested a vdW functional.")

    return calc
