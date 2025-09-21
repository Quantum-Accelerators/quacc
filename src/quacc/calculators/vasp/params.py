"""Parameter-related utilities for the Vasp calculator."""

from __future__ import annotations

from importlib.util import find_spec
from logging import getLogger
from typing import TYPE_CHECKING

import numpy as np
import psutil
from ase.calculators.vasp import Vasp as Vasp_
from monty.dev import requires
from pymatgen.io.ase import AseAtomsAdaptor

from quacc.atoms.core import check_is_metal
from quacc.utils.dicts import sort_dict
from quacc.utils.kpts import convert_pmg_kpts

has_atomate2 = bool(find_spec("atomate2"))

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase.atoms import Atoms
    from pymatgen.io.vasp.sets import DictSet

    from quacc.types import PmgKpts, SourceDirectory

    if has_atomate2:
        from atomate2.vasp.jobs.base import BaseVaspMaker

LOGGER = getLogger(__name__)


def get_param_swaps(
    user_calc_params: dict[str, Any],
    pmg_kpts: dict[Literal["line_density", "kppvol", "kppa"], float],
    input_atoms: Atoms,
    incar_copilot: Literal["off", "on", "aggressive"],
) -> dict[str, Any]:
    """
    Swaps out bad INCAR flags.

    Parameters
    ----------
    user_calc_params
        The user-provided calculator parameters.
    pmg_kpts
        The pmg_kpts kwarg.
    input_atoms
        The input atoms.
    incar_copilot
        INCAR copilot mode. See `quacc.calculators.vasp.vasp.Vasp` for more info.

    Returns
    -------
    dict
        The updated user-provided calculator parameters.
    """
    is_metal = check_is_metal(input_atoms)
    calc = Vasp_(**user_calc_params)
    max_Z = input_atoms.get_atomic_numbers().max()

    if calc.parameters.get("lmaxmix", 2) < 6 and max_Z > 56:
        LOGGER.info("Recommending LMAXMIX = 6 because you have f electrons.")
        calc.set(lmaxmix=6)
    elif calc.parameters.get("lmaxmix", 2) < 4 and max_Z > 20:
        LOGGER.info("Recommending LMAXMIX = 4 because you have d electrons.")
        calc.set(lmaxmix=4)

    if (
        calc.parameters.get("luse_vdw", False)
        or calc.parameters.get("lhfcalc", False)
        or calc.parameters.get("ldau", False)
        or calc.parameters.get("ldau_luj", {})
        or calc.parameters.get("metagga", "")
    ) and not calc.parameters.get("lasph"):
        LOGGER.info(
            "Recommending LASPH = True because you have a +U, vdW, meta-GGA, or hybrid calculation."
        )
        calc.set(lasph=True)

    if calc.parameters.get("metagga", "") and (
        calc.parameters.get("algo", "normal").lower() != "all"
    ):
        LOGGER.info("Recommending ALGO = All because you have a meta-GGA calculation.")
        calc.set(algo="all", isearch=1)

    if calc.parameters.get("lhfcalc", False) and (
        calc.parameters.get("algo", "normal").lower() != "normal"
    ):
        LOGGER.info("Recommending ALGO = Normal because you have a hybrid calculation.")
        calc.set(algo="normal")

    if (
        is_metal
        and (calc.parameters.get("ismear", 1) < 0)
        and (calc.parameters.get("nsw", 0) > 0)
    ):
        LOGGER.info(
            "Recommending ISMEAR = 1 and SIGMA = 0.1 because you are likely relaxing a metal."
        )
        calc.set(ismear=1, sigma=0.1)

    if (
        calc.parameters.get("ismear", 1) != -5
        and calc.parameters.get("nsw", 0) == 0
        and (
            (calc.kpts is not None and np.prod(calc.kpts) >= 4)
            or (calc.float_params["kspacing"] and calc.float_params["kspacing"] <= 0.5)
        )
    ):
        LOGGER.info("Recommending ISMEAR = -5 because you have a static calculation.")
        calc.set(ismear=-5)

    if (
        calc.parameters.get("ismear", 1) == -5
        and (calc.kpts is not None and np.prod(calc.kpts) < 4)
        and calc.parameters.get("kspacing", None) is None
    ):
        LOGGER.info(
            "Recommending ISMEAR = 0 because you don't have enough k-points for ISMEAR = -5."
        )
        calc.set(ismear=0)

    if (
        calc.parameters.get("kspacing", 0.5) > 0.5
        and calc.parameters.get("ismear", 1) == -5
    ):
        LOGGER.info(
            "Recommending ISMEAR = 0 because KSPACING is likely too large for ISMEAR = -5."
        )
        calc.set(ismear=0)

    if (
        pmg_kpts
        and pmg_kpts.get("line_density")
        and calc.parameters.get("ismear", 1) != 0
    ):
        LOGGER.info(
            "Recommending ISMEAR = 0 and SIGMA = 0.01 because you are doing a line mode calculation."
        )
        calc.set(ismear=0, sigma=0.01)

    if calc.parameters.get("ismear", 1) == 0 and (
        calc.parameters.get("sigma", 0.2) > 0.05
    ):
        LOGGER.info(
            "Recommending SIGMA = 0.05 because ISMEAR = 0 was requested with SIGMA > 0.05."
        )
        calc.set(sigma=0.05)

    if calc.parameters.get("nsw", 0) > 0 and calc.parameters.get("laechg", False):
        LOGGER.info(
            "Recommending LAECHG = False because you have NSW > 0. LAECHG is not compatible with NSW > 0."
        )
        calc.set(laechg=False)

    if calc.parameters.get("ldauprint", 0) == 0 and (
        calc.parameters.get("ldau", False) or calc.parameters.get("ldau_luj", {})
    ):
        LOGGER.info("Recommending LDAUPRINT = 1 because LDAU = True.")
        calc.set(ldauprint=1)

    if calc.parameters.get("lreal", False) and len(input_atoms) < 30:
        LOGGER.info(
            "Recommending LREAL = False because you have a small system (< 30 atoms/cell)."
        )
        calc.set(lreal=False)

    if not calc.parameters.get("lorbit", False) and (
        calc.parameters.get("ispin", 1) == 2
        or np.any(input_atoms.get_initial_magnetic_moments() != 0)
    ):
        LOGGER.info(
            "Recommending LORBIT = 11 because you have a spin-polarized calculation."
        )
        calc.set(lorbit=11)

    if not calc.parameters.get("npar") and not calc.parameters.get("ncore"):
        ncores = psutil.cpu_count(logical=False) or 1
        for ncore in range(int(np.sqrt(ncores)), ncores):
            if ncores % ncore == 0:
                LOGGER.info(
                    f"Recommending NCORE = {ncore} per the sqrt(# cores) suggestion by VASP."
                )
                calc.set(ncore=ncore, npar=None)
                break

    if (
        calc.parameters.get("ncore", 1) > 1
        or (calc.parameters.get("npar") and calc.parameters.get("npar", 1) > 1)
    ) and (
        calc.parameters.get("lhfcalc", False) is True
        or calc.parameters.get("lrpa", False) is True
        or calc.parameters.get("lepsilon", False) is True
        or calc.parameters.get("ibrion", 0) in [5, 6, 7, 8]
    ):
        LOGGER.info(
            "Recommending NCORE = 1 because NCORE/NPAR is not compatible with this job type."
        )
        calc.set(ncore=1, npar=None)

    if (
        calc.parameters.get("kpar")
        and (
            calc.kpts is not None
            and calc.parameters.get("kpar", 1) > np.prod(calc.kpts)
        )
        and calc.float_params["kspacing"] is None
    ):
        LOGGER.info(
            "Recommending KPAR = 1 because you have too few k-points to parallelize."
        )
        calc.set(kpar=1)

    if (
        calc.parameters.get("lhfcalc", False) is True
        and calc.parameters.get("isym", 3) < 3
    ):
        LOGGER.info(
            "Recommending ISYM = 3 because you are running a hybrid calculation."
        )
        calc.set(isym=3)

    if calc.parameters.get("lsorbit", False):
        LOGGER.info(
            "Recommending ISYM = -1 because you are running an SOC calculation."
        )
        calc.set(isym=-1)

    if (
        calc.parameters.get("algo", "normal") in ("all", "conjugate")
        and calc.parameters.get("isearch", 0) != 1
    ):
        LOGGER.info("Recommending ISEARCH = 1 because you have ALGO = All.")
        calc.set(isearch=1)

    if (
        calc.parameters.get("isif", 2) in (3, 6, 7, 8)
        and calc.parameters.get("nsw", 0) > 0
    ):
        if calc.encut is None:
            LOGGER.warning(
                "Be careful of Pulay stresses. At the end of your run, re-relax your structure with your current ENCUT or set ENCUT=1.3*max(ENMAX)."
            )
        if "He" in input_atoms.get_chemical_symbols() and (
            calc.encut is None or calc.encut < 478.896 * 1.3
        ):
            LOGGER.warning(
                "Be careful of Pulay stresses. At the end of your run, re-relax your structure with your current ENCUT or set ENCUT>=623."
            )

        if (
            "Li" in input_atoms.get_chemical_symbols()
            and calc.parameters.get("setups", {})
            and isinstance(calc.parameters["setups"], dict)
            and calc.parameters["setups"].get("Li", "") in ("Li_sv", "_sv")
            and (calc.encut is None or calc.encut < 499.034 * 1.3)
        ):
            LOGGER.warning(
                "Be careful of Pulay stresses. At the end of your run, re-relax your structure with your current ENCUT or set ENCUT>=650."
            )

    if (
        calc.parameters.get("metagga", "").lower() == "r2scan"
        and calc.parameters.get("ivdw", 0) == 13
        and not calc.parameters.get("vdw_s6")
        and not calc.parameters.get("vdw_s8")
        and not calc.parameters.get("vdw_a1")
        and not calc.parameters.get("vdw_a2")
    ):
        LOGGER.info("Setting VDW_S6, VDW_S8, VDW_A1, VDW_A2 parameters for r2SCAN-D4.")
        calc.set(vdw_s6=1.0, vdw_s8=0.60187490, vdw_a1=0.51559235, vdw_a2=5.77342911)

    if (
        calc.parameters.get("lhfcalc", False)
        and calc.parameters.get("hfscreen", 0) == 0.2
        and calc.parameters.get("ivdw", 1) == 12
        and not calc.parameters.get("vdw_s6")
        and not calc.parameters.get("vdw_s8")
        and not calc.parameters.get("vdw_a1")
        and not calc.parameters.get("vdw_a2")
    ):
        LOGGER.info(
            "Setting VDW_S6, VDW_S8, VDW_A1, VDW_A2 parameters for HSE06-D3(BJ)."
        )
        calc.set(vdw_s6=1.0, vdw_s8=2.310, vdw_a1=0.383, vdw_a2=5.685)

    if (
        input_atoms.get_chemical_formula() == "O2"
        and input_atoms.get_initial_magnetic_moments().sum() == 0
    ):
        LOGGER.warning(
            "You are running O2 without magnetic moments, but its ground state should have 2 unpaired electrons!"
        )

    new_parameters = (
        calc.parameters
        if incar_copilot == "aggressive"
        else (
            calc.parameters | user_calc_params
            if incar_copilot == "on"
            else user_calc_params
        )
    )
    if changed_parameters := {
        k: new_parameters[k] for k in set(new_parameters) - set(user_calc_params)
    }:
        LOGGER.info(
            f"The following parameters were changed: {sort_dict(changed_parameters)}"
        )

    return new_parameters


def remove_unused_flags(user_calc_params: dict[str, Any]) -> dict[str, Any]:
    """
    Removes unused flags in the INCAR, like EDIFFG if you are doing NSW = 0.

    Parameters
    ----------
    user_calc_params
        The updated user-provided calculator parameters.

    Returns
    -------
    dict
        The updated user-provided calculator parameters.
    """
    if user_calc_params.get("nsw", 0) == 0:
        # Turn off opt flags if NSW = 0
        opt_flags = ("ediffg", "ibrion", "potim", "iopt")
        for opt_flag in opt_flags:
            user_calc_params.pop(opt_flag, None)

    if not user_calc_params.get("ldau", False) and not user_calc_params.get("ldau_luj"):
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
        for ldau_flag in ldau_flags:
            user_calc_params.pop(ldau_flag, None)

    # Handle kspacing flags
    if user_calc_params.get("kspacing"):
        user_calc_params["gamma"] = None
        user_calc_params["kpts"] = None
    else:
        user_calc_params.pop("kgamma", None)

    # Remove None keys
    none_keys = [
        k
        for k, v in user_calc_params.items()
        if v is None and k not in Vasp_().input_params
    ]
    for none_key in none_keys:
        del user_calc_params[none_key]

    return user_calc_params


def normalize_params(user_calc_params: dict[str, Any]) -> dict[str, Any]:
    """
    Normalizes the user-provided calculator parameters.

    Parameters
    ----------
    user_calc_params
        The user-provided calculator parameters.

    Returns
    -------
    dict
        The updated user-provided calculator parameters.
    """
    for k, v in user_calc_params.items():
        if isinstance(v, str):
            user_calc_params[k] = v.lower()
    return user_calc_params


def set_auto_dipole(
    user_calc_params: dict[str, Any], input_atoms: Atoms
) -> dict[str, Any]:
    """
    Sets flags related to the auto_dipole kwarg.

    Parameters
    ----------
    user_calc_params
        The user-provided calculator parameters.
    input_atoms
        The input atoms.

    Returns
    -------
    dict
        The updated user-provided calculator parameters.
    """
    com = input_atoms.get_center_of_mass(scaled=True)
    if "dipol" not in user_calc_params:
        user_calc_params["dipol"] = com
    if "idipol" not in user_calc_params:
        user_calc_params["idipol"] = 3
    if "ldipol" not in user_calc_params:
        user_calc_params["ldipol"] = True

    return user_calc_params


def set_pmg_kpts(
    user_calc_params: PmgKpts,
    pmg_kpts: dict[Literal["line_density", "kppvol", "kppa"], float],
    input_atoms: Atoms,
) -> dict[str, Any]:
    """
    Shortcuts for pymatgen k-point generation schemes.

    Parameters
    ----------
    user_calc_params
        The user-provided calculator parameters.
    pmg_kpts
        The pmg_kpts kwarg.
    input_atoms
        The input atoms.

    Returns
    -------
    dict
        The updated user-provided calculator parameters.
    """
    kpts, gamma = convert_pmg_kpts(
        pmg_kpts, input_atoms, force_gamma=user_calc_params.get("gamma", False)
    )
    reciprocal = bool(pmg_kpts.get("line_density"))

    user_calc_params["kpts"] = kpts
    if reciprocal and user_calc_params.get("reciprocal") is None:
        user_calc_params["reciprocal"] = reciprocal
    if user_calc_params.get("gamma") is None:
        user_calc_params["gamma"] = gamma

    return user_calc_params


class MPtoASEConverter:
    """
    Convert an MP-formatted input set to an ASE-formatted input set.
    """

    def __init__(
        self, atoms: Atoms | None = None, prev_dir: SourceDirectory | None = None
    ) -> None:
        """
        Initialize the converter.

        Parameters
        ----------
        atoms
            The ASE atoms object.
        prev_dir
            The previous directory.

        Returns
        -------
        None
        """
        if atoms is None and prev_dir is None:
            raise ValueError("Either atoms or prev_dir must be provided.")
        self.atoms = atoms
        self.prev_dir = prev_dir
        if self.atoms:
            self.ase_sort, self.ase_resort = Vasp_()._make_sort(self.atoms)
            self.structure = AseAtomsAdaptor.get_structure(self.atoms[self.ase_sort])
        else:
            self.structure = None

    def convert_dict_set(self, dict_set: DictSet) -> dict:
        """
        Convert a Pymatgen DictSet to a dictionary of ASE VASP parameters.

        Parameters
        ----------
        dict_set
            The instantiated Pymatgen DictSet.

        Returns
        -------
        dict
            The ASE VASP parameters.
        """
        assert hasattr(dict_set, "sort_structure")
        dict_set.sort_structure = False
        vasp_input = dict_set.get_input_set(
            structure=self.structure, potcar_spec=True, prev_dir=self.prev_dir
        )
        self.incar_dict = vasp_input["INCAR"]
        self.pmg_kpts = vasp_input.get("KPOINTS")
        self.potcar_symbols = vasp_input["POTCAR.spec"].split("\n")
        self.potcar_functional = dict_set.potcar_functional
        self.poscar = vasp_input["POSCAR"]
        return self._convert()

    @requires(has_atomate2, "atomate2 is not installed.")
    def convert_vasp_maker(self, VaspMaker: BaseVaspMaker) -> dict:
        """
        Convert an atomate2 VaspMaker to a dictionary of ASE VASP parameters.

        Parameters
        ----------
        VaspMaker
            The instantiated atomate2 VaspMaker.

        Returns
        -------
        dict
            The ASE VASP parameters.
        """
        input_set_generator = VaspMaker.input_set_generator
        assert hasattr(input_set_generator, "sort_structure")
        input_set_generator.sort_structure = False
        input_set = input_set_generator.get_input_set(
            structure=self.structure, potcar_spec=True, prev_dir=self.prev_dir
        )
        self.incar_dict = input_set.incar
        self.pmg_kpts = input_set.kpoints
        self.potcar_symbols = (
            input_set.potcar.split("\n")
            if isinstance(input_set.potcar, str)
            else input_set.potcar
        )
        self.potcar_functional = input_set_generator.potcar_functional
        self.poscar = input_set.poscar
        return self._convert()

    def _convert(self) -> dict:
        """
        Convert the MP input to a dictionary of ASE VASP parameters.

        Returns
        -------
        dict
            The ASE VASP parameters.
        """
        self.incar_dict = {k.lower(): v for k, v in self.incar_dict.items()}
        pp = self.potcar_functional.split("_")[0]
        assert pp.lower() in ["lda", "pw91", "pbe"]
        potcar_setups = {symbol.split("_")[0]: symbol for symbol in self.potcar_symbols}
        for k, v in potcar_setups.items():
            if k in v:
                potcar_setups[k] = v.split(k)[-1]

        full_input_params = self.incar_dict | {"setups": potcar_setups, "pp": pp}

        if self.pmg_kpts:
            kpts_dict = self.pmg_kpts.as_dict()
            full_input_params |= {
                "kpts": kpts_dict["kpoints"][0],
                "gamma": kpts_dict["generation_style"].lower() == "gamma",
            }

        return full_input_params
