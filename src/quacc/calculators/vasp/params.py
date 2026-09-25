"""Parameter-related utilities for the Vasp calculator."""

from __future__ import annotations

from dataclasses import dataclass, field
from importlib.util import find_spec
from logging import INFO, WARNING, getLogger
from typing import TYPE_CHECKING

import numpy as np
import psutil
from ase.calculators.vasp import Vasp as Vasp_
from monty.dev import requires
from pymatgen.io.ase import AseAtomsAdaptor

from quacc.atoms.core import check_is_metal
from quacc.utils.kpts import convert_pmg_kpts

has_atomate2 = bool(find_spec("atomate2"))

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase.atoms import Atoms
    from pymatgen.io.vasp.sets import VaspInputSet

    from quacc.types import PmgKpts, SourceDirectory

    if has_atomate2:
        from atomate2.vasp.jobs.base import BaseVaspMaker
        from atomate2.vasp.sets.base import VaspInputGenerator

LOGGER = getLogger(__name__)


@dataclass(frozen=True)
class CopilotChange:
    """An INCAR parameter change made by the copilot."""

    parameter: str
    old_value: Any
    new_value: Any
    reason: str


@dataclass(frozen=True)
class CopilotWarning:
    """An issue that the INCAR copilot cannot safely fix."""

    message: str
    remediation: str | None = None
    reason: str | None = None


@dataclass
class CopilotReport:
    """Structured results from the INCAR copilot."""

    changes: list[CopilotChange] = field(default_factory=list)
    warnings: list[CopilotWarning] = field(default_factory=list)


def log_copilot_report(report: CopilotReport) -> None:
    """Log a framed summary as a single record."""
    header = "================ VASP INCAR COPILOT ================"
    footer = "=" * len(header)
    lines = [header]
    if report.warnings:
        lines.append("⚠️ Attention required:")
        for warning in report.warnings:
            lines.append(f"  {warning.message}")
            if warning.reason:
                lines.append(f"    Reason: {warning.reason}")
            if warning.remediation:
                lines.append(f"    Recommendation: {warning.remediation}")
    if report.changes:
        if report.warnings:
            lines.append("")
        lines.append("Applied changes:")
        lines.extend(
            f"  {change.parameter.upper()}: {change.old_value!r} -> {change.new_value!r}\n"
            f"    Reason: {change.reason}"
            for change in report.changes
        )
    lines.append(footer)

    log_level = WARNING if report.warnings else INFO
    LOGGER.log(log_level, "\n%s", "\n".join(lines))


def get_param_swaps(
    user_calc_params: dict[str, Any],
    input_atoms: Atoms,
    pmg_kpts: dict[Literal["line_density", "kppvol", "kppa"], float] | None = None,
    incar_copilot_mode: Literal[
        "off", "critical", "standard", "aggressive"
    ] = "standard",
) -> dict[str, Any]:
    """
    Swaps out bad INCAR flags.

    Parameters
    ----------
    user_calc_params
        The user-provided calculator parameters.
    input_atoms
        The input atoms.
    pmg_kpts
        The pmg_kpts kwarg.
    incar_copilot_mode
        INCAR copilot mode. See `quacc.calculators.vasp.vasp.Vasp` for more info.

    Returns
    -------
    dict
        The updated user-provided calculator parameters.
    """
    report = CopilotReport()
    is_metal = check_is_metal(input_atoms)
    calc = Vasp_(**remove_unused_flags(user_calc_params))
    change_reasons: dict[str, str] = {}

    def recommend(reason: str, **parameters: Any) -> None:
        change_reasons.update(dict.fromkeys(parameters, reason))
        calc.set(**parameters)

    max_Z = input_atoms.get_atomic_numbers().max()

    # ----------------------------
    # General INCAR swaps
    # ----------------------------
    if incar_copilot_mode.lower() not in {"off", "critical"}:
        if calc.parameters.get("lmaxmix", 2) < 6 and max_Z > 56:
            recommend("f electrons were detected.", lmaxmix=6)
        elif calc.parameters.get("lmaxmix", 2) < 4 and max_Z > 20:
            recommend("d electrons were detected.", lmaxmix=4)

        if (
            calc.parameters.get("luse_vdw", False)
            or calc.parameters.get("lhfcalc", False)
            or calc.parameters.get("ldau", False)
            or calc.parameters.get("ldau_luj", {})
            or calc.parameters.get("metagga", "")
        ) and not calc.parameters.get("lasph"):
            recommend(
                "LASPH is recommended for +U, vdW, meta-GGA, and hybrid calculations.",
                lasph=True,
            )

        if calc.parameters.get("metagga", "") and (
            calc.parameters.get("algo", "normal").lower() != "all"
        ):
            recommend(
                "ALGO=All is recommended for meta-GGA calculations.",
                algo="all",
                isearch=1,
            )

        if (
            is_metal
            and (calc.parameters.get("ismear", 1) < 0)
            and (calc.parameters.get("nsw", 0) > 0)
        ):
            recommend(
                "Metal relaxations should use finite-width smearing.",
                ismear=1,
                sigma=0.1,
            )

        if (
            pmg_kpts
            and pmg_kpts.get("line_density")
            and calc.parameters.get("ismear", 1) != 0
        ):
            recommend(
                "Line-mode calculations should use Gaussian smearing.",
                ismear=0,
                sigma=0.01,
            )

        if calc.parameters.get("ismear", 1) == 0 and (
            calc.parameters.get("sigma", 0.2) > 0.05
        ):
            recommend("SIGMA should be <=0.05 when ISMEAR=0.", sigma=0.05)

        if calc.parameters.get("nsw", 0) > 0 and calc.parameters.get("laechg", False):
            recommend("LAECHG is not compatible with NSW>0.", laechg=False)

        if calc.parameters.get("ldauprint", 0) == 0 and (
            calc.parameters.get("ldau", False) or calc.parameters.get("ldau_luj", {})
        ):
            recommend("LDAUPRINT=1 is recommended when LDAU is enabled.", ldauprint=1)

        if calc.parameters.get("lreal", False) and len(input_atoms) < 30:
            recommend(
                "Reciprocal-space projectors are recommended for systems with fewer than 30 atoms.",
                lreal=False,
            )

        if (
            calc.parameters.get("lhfcalc", False) is True
            and calc.parameters.get("isym", 3) < 3
        ):
            recommend("ISYM=3 is recommended for hybrid calculations.", isym=3)

        if calc.parameters.get("lsorbit", False):
            recommend("Symmetry should be disabled for spin-orbit coupling.", isym=-1)

        if (
            calc.parameters.get("algo", "normal") in ("all", "conjugate")
            and calc.parameters.get("isearch", 0) != 1
        ):
            recommend("ISEARCH=1 is required when ALGO=All.", isearch=1)

    # ----------------------------
    # Critical INCAR swaps
    # ----------------------------
    pre_critical_params = dict(calc.parameters)
    if incar_copilot_mode.lower() != "off":
        if (ediff := calc.parameters.get("ediff", 1e-4)) > 1e-4:
            report.warnings.append(
                CopilotWarning(
                    f"EDIFF={ediff!r} is greater than 1e-4; the results are likely unconverged."
                )
            )

        if not calc.parameters.get("lorbit", False) and (
            calc.parameters.get("ispin", 1) == 2
            or np.any(input_atoms.get_initial_magnetic_moments() != 0)
        ):
            recommend(
                "LORBIT=11 is recommended for spin-polarized calculations.", lorbit=11
            )

        if (
            calc.parameters.get("ismear", 1) == -5
            and (calc.kpts is not None and np.prod(calc.kpts) < 4)
            and calc.parameters.get("kspacing", None) is None
        ):
            recommend("There are too few k-points for tetrahedron smearing.", ismear=0)

        if (
            calc.parameters.get("kspacing", 0.5) > 0.5
            and calc.parameters.get("ismear", 1) == -5
        ):
            recommend("KSPACING is too large for tetrahedron smearing.", ismear=0)

        if calc.parameters.get("ismear", 1) == -5 and calc.parameters.get(
            "algo", "normal"
        ) in ("all", "conjugate", "damped"):
            recommend(
                "The selected ALGO is incompatible with ISMEAR=-5.", algo="normal"
            )

        if calc.parameters.get("lhfcalc", False) and (
            calc.parameters.get("algo", "normal").lower() != "normal"
        ):
            recommend("ALGO=Normal is required for hybrid calculations.", algo="normal")

        if (
            calc.parameters.get("ncore", 1) > 1
            or (calc.parameters.get("npar") and calc.parameters.get("npar", 1) > 1)
        ) and (
            calc.parameters.get("lhfcalc", False) is True
            or calc.parameters.get("lrpa", False) is True
            or calc.parameters.get("lepsilon", False) is True
            or calc.parameters.get("ibrion", 0) in [5, 6, 7, 8]
        ):
            recommend(
                "NCORE/NPAR is incompatible with this job type.", ncore=1, npar=None
            )

        if not calc.parameters.get("npar") and not calc.parameters.get("ncore"):
            ncores = psutil.cpu_count(logical=False) or 1
            for ncore in range(int(np.sqrt(ncores)), ncores):
                if ncores % ncore == 0:
                    recommend(
                        "This follows VASP's sqrt(number of cores) recommendation.",
                        ncore=ncore,
                        npar=None,
                    )
                    break

        if (
            calc.parameters.get("kpar")
            and (
                calc.kpts is not None
                and calc.parameters.get("kpar", 1) > np.prod(calc.kpts)
            )
            and calc.float_params["kspacing"] is None
        ):
            recommend(
                "There are too few k-points to parallelize over k-points.", kpar=1
            )

        if (
            calc.parameters.get("isif", 2) in (3, 6, 7, 8)
            and calc.parameters.get("nsw", 0) > 0
        ):
            if calc.encut is None:
                report.warnings.append(
                    CopilotWarning(
                        "Pulay stress risk during variable-cell relaxation because ENCUT was not explicitly set.",
                        "Re-relax the final structure with the converged ENCUT or set ENCUT=1.3*max(ENMAX).",
                    )
                )
            if "He" in input_atoms.get_chemical_symbols() and (
                calc.encut is None or calc.encut < 478.896 * 1.3
            ):
                report.warnings.append(
                    CopilotWarning(
                        "Pulay stress risk for a variable-cell relaxation containing He.",
                        "Re-relax the final structure with the converged ENCUT or set ENCUT>=623.",
                    )
                )

            if (
                "Li" in input_atoms.get_chemical_symbols()
                and calc.parameters.get("setups", {})
                and isinstance(calc.parameters["setups"], dict)
                and calc.parameters["setups"].get("Li", "") in ("Li_sv", "_sv")
                and (calc.encut is None or calc.encut < 499.034 * 1.3)
            ):
                report.warnings.append(
                    CopilotWarning(
                        "Pulay stress risk for a variable-cell relaxation using Li_sv.",
                        "Re-relax the final structure with the converged ENCUT or set ENCUT>=650.",
                    )
                )

        if (
            (
                calc.parameters.get("xc", "").lower() == "r2scan"
                or calc.parameters.get("metagga", "").lower() == "r2scan"
            )
            and calc.parameters.get("ivdw", 0) == 13
            and not calc.parameters.get("vdw_s6")
            and not calc.parameters.get("vdw_s8")
            and not calc.parameters.get("vdw_a1")
            and not calc.parameters.get("vdw_a2")
        ):
            recommend(
                "These are the recommended r2SCAN-D4 damping parameters.",
                vdw_s6=1.0,
                vdw_s8=0.60187490,
                vdw_a1=0.51559235,
                vdw_a2=5.77342911,
            )

        if (
            calc.parameters.get("lhfcalc", False)
            and calc.parameters.get("hfscreen", 0) == 0.2
            and calc.parameters.get("ivdw", 1) == 12
            and not calc.parameters.get("vdw_s6")
            and not calc.parameters.get("vdw_s8")
            and not calc.parameters.get("vdw_a1")
            and not calc.parameters.get("vdw_a2")
        ):
            recommend(
                "These are the recommended HSE06-D3(BJ) damping parameters.",
                vdw_s6=1.0,
                vdw_s8=2.310,
                vdw_a1=0.383,
                vdw_a2=5.685,
            )

        if input_atoms.get_chemical_formula() == "O2" and all(
            input_atoms.get_initial_magnetic_moments() == 0
        ):
            report.warnings.append(
                CopilotWarning(
                    "O2 is being run without magnetic moments, but its ground state should have 2 unpaired electrons."
                )
            )

    # ----------------------------
    # Finalize INCAR swaps
    # ----------------------------
    critical_swap_changes = {
        k: v
        for k, v in calc.parameters.items()
        if _params_differ(pre_critical_params.get(k), v)
    }
    recommended_params = dict(calc.parameters)

    if incar_copilot_mode == "aggressive":
        new_parameters = calc.parameters
    else:
        new_parameters = (calc.parameters | user_calc_params) | critical_swap_changes

    report.changes = [
        CopilotChange(
            k,
            user_calc_params.get(k),
            new_parameters[k],
            change_reasons.get(k, "General INCAR compatibility rule was applied."),
        )
        for k in sorted(new_parameters)
        if k in change_reasons
        and _params_differ(user_calc_params.get(k), new_parameters[k])
    ]

    overridden_user_params = {
        k: (user_calc_params[k], new_parameters[k])
        for k in user_calc_params
        if k in new_parameters
        and _params_differ(new_parameters[k], user_calc_params[k])
    }
    for k, (old, new) in overridden_user_params.items():
        report.warnings.append(
            CopilotWarning(
                f"{k.upper()} was changed from {old!r} to {new!r}.",
                reason=change_reasons.get(
                    k, "General INCAR compatibility rule was applied."
                ),
            )
        )

    if overridden_swaps := {
        k: (user_calc_params.get(k), recommended_params[k])
        for k in recommended_params
        if _params_differ(recommended_params[k], new_parameters.get(k))
    }:
        for k, (current, recommended) in overridden_swaps.items():
            report.warnings.append(
                CopilotWarning(
                    f"{k.upper()} was not changed from {current!r} to {recommended!r}, but should be modified.",
                    reason=change_reasons.get(k),
                )
            )

    if incar_copilot_mode.lower() != "off":
        log_copilot_report(report)

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
    new_user_calc_params = user_calc_params.copy()
    if new_user_calc_params.get("nsw", 0) == 0:
        # Turn off opt flags if NSW = 0
        opt_flags = ("ediffg", "ibrion", "potim", "iopt")
        for opt_flag in opt_flags:
            new_user_calc_params.pop(opt_flag, None)

    if not new_user_calc_params.get("ldau", False) and not new_user_calc_params.get(
        "ldau_luj"
    ):
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
            new_user_calc_params.pop(ldau_flag, None)

    # Handle kspacing flags
    if new_user_calc_params.get("kspacing"):
        new_user_calc_params["gamma"] = None
        new_user_calc_params["kpts"] = None
    else:
        new_user_calc_params.pop("kgamma", None)

    # Remove None keys
    none_keys = [
        k
        for k, v in new_user_calc_params.items()
        if v is None and k not in Vasp_().input_params
    ]
    for none_key in none_keys:
        del new_user_calc_params[none_key]

    return new_user_calc_params


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
    return {
        (k.lower() if isinstance(k, str) else k): (
            v.lower() if isinstance(v, str) else v
        )
        for k, v in user_calc_params.items()
    }


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

    def convert_input_set(self, input_set: VaspInputSet) -> dict:
        """
        Convert a Pymatgen VaspInputSet to a dictionary of ASE VASP parameters.

        Parameters
        ----------
        input_set
            The instantiated Pymatgen VaspInputSet.

        Returns
        -------
        dict
            The ASE VASP parameters.
        """
        assert hasattr(input_set, "sort_structure")
        input_set.sort_structure = False
        vasp_input = input_set.get_input_set(
            structure=self.structure, potcar_spec=True, prev_dir=self.prev_dir
        )
        self.incar_dict = vasp_input["INCAR"]
        self.pmg_kpts = vasp_input.get("KPOINTS")
        self.potcar_symbols = vasp_input["POTCAR.spec"].split("\n")
        self.potcar_functional = input_set.potcar_functional
        self.poscar = vasp_input["POSCAR"]
        return self._convert()

    @requires(has_atomate2, "atomate2 is not installed.")
    def convert_input_generator(self, input_generator: VaspInputGenerator) -> dict:
        """
        Convert a VaspInputGenerator to a dictionary of ASE VASP parameters.

        Parameters
        ----------
        input_generator
            The instantiated VaspInputGenerator.

        Returns
        -------
        dict
            The ASE VASP parameters.
        """
        assert hasattr(input_generator, "sort_structure")
        input_generator.sort_structure = False
        input_set = input_generator.get_input_set(
            structure=self.structure, potcar_spec=True, prev_dir=self.prev_dir
        )
        self.incar_dict = input_set.incar
        self.pmg_kpts = input_set.kpoints
        self.potcar_symbols = (
            input_set.potcar.split("\n")
            if isinstance(input_set.potcar, str)
            else input_set.potcar
        )
        self.potcar_functional = input_generator.potcar_functional
        self.poscar = input_set.poscar
        return self._convert()

    @requires(has_atomate2, "atomate2 is not installed.")
    def convert_maker(self, VaspMaker: BaseVaspMaker) -> dict:
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
        return self.convert_input_generator(VaspMaker.input_set_generator)

    def _convert(self) -> dict:
        """
        Convert the MP input to a dictionary of ASE VASP parameters.

        Returns
        -------
        dict
            The ASE VASP parameters.
        """
        self.incar_dict = {k.lower(): v for k, v in self.incar_dict.items()}
        parts = self.potcar_functional.split("_")
        pp = parts[0]
        pp_version = parts[1] if len(parts) > 1 else "original"
        assert pp.lower() in ["lda", "pbe"]
        potcar_setups = {symbol.split("_")[0]: symbol for symbol in self.potcar_symbols}
        for k, v in potcar_setups.items():
            if k in v:
                potcar_setups[k] = v.split(k)[-1]

        full_input_params = self.incar_dict | {
            "setups": potcar_setups,
            "pp": pp,
            "pp_version": pp_version,
        }

        if self.pmg_kpts:
            kpts_dict = self.pmg_kpts.as_dict()
            full_input_params |= {
                "kpts": kpts_dict["kpoints"][0],
                "gamma": kpts_dict["generation_style"].lower() == "gamma",
            }

        return full_input_params


def _params_differ(a: object, b: object) -> bool:
    """Compare two parameter values, handling NumPy arrays."""
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        return not np.array_equal(a, b)
    elif isinstance(a, str) and isinstance(b, str):
        return a.lower() != b.lower()
    return a != b
