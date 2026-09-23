"""Meta FAIR recipes"""

from __future__ import annotations

from importlib.util import find_spec
from math import ceil
from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import job
from quacc.calculators.vasp.params import MPtoASEConverter
from quacc.recipes.vasp._base import run_and_summarize

has_fairchem = bool(find_spec("fairchem"))
has_fairchem_omat = (
    has_fairchem
    and bool(find_spec("fairchem.data"))
    and bool(find_spec("fairchem.data.omat"))
)
has_fairchem_oc = (
    has_fairchem
    and bool(find_spec("fairchem.data"))
    and bool(find_spec("fairchem.data.oc"))
)
has_atomate2 = bool(find_spec("atomate2"))

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import CopyFiles, VaspSchema


@job
@requires(
    has_fairchem_omat,
    "fairchem-data-omat is not installed. Run `pip install quacc[fairchem]`",
)
def omat_static_job(
    atoms: Atoms,
    copy_files: CopyFiles | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs: Any,
) -> VaspSchema:
    """
    Carry out a static calculation with OMat settings.

    Parameters
    ----------
    atoms
        Atoms object
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """
    from fairchem.data.omat.vasp.sets import OMat24StaticSet

    calc_defaults = MPtoASEConverter(atoms=atoms).convert_input_set(OMat24StaticSet())
    calc_defaults |= {"pp_version": "54", "incar_copilot_mode": "critical"}

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "OMat Static"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[fairchem]`")
def omc_static_job(
    atoms: Atoms,
    copy_files: CopyFiles | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs: Any,
) -> VaspSchema:
    """
    Carry out a static calculation with OMC settings.

    Parameters
    ----------
    atoms
        Atoms object
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """

    calc_defaults = _make_omc_inputs(atoms)
    calc_defaults |= {
        "pp_version": "54",
        "incar_copilot_mode": "critical",
        "use_custodian": False,
    }

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "OMC Static"} | (additional_fields or {}),
        copy_files=copy_files,
    )


def _make_omc_inputs(atoms: Atoms) -> dict:
    """
    Helper function to make a fairchem input set.

    Parameters
    ----------
    atoms
        Atoms object
    dataset
        Dataset to use. Currently only "omc" is supported.

    Returns
    -------
    dict
        Dictionary of ASE VASP calculator parameters.
    """
    from atomate2.vasp.sets.core import StaticSetGenerator

    input_generator = StaticSetGenerator(
        user_incar_settings={
            "ADDGRID": True,
            "ALGO": "Normal",
            "EDIFF": 1e-06,
            "ENCUT": 520,
            "GGA": "PE",
            "IBRION": -1,
            "ISIF": 0,
            "ISMEAR": 0,
            "ISPIN": 1,
            "IVDW": 11,
            "LREAL": False,
            "LMIXTAU": True,
            "LASPH": True,
            "LORBIT": 11,
            "LWAVE": False,
            "LAECHG": False,
            "LVTOT": False,
            "NELM": 200,
            "NELMDL": -10,
            "NSW": 0,
            "PREC": "Normal",
            "SIGMA": 0.1,
        },
        user_potcar_functional="PBE_54_W_HASH",
        auto_kspacing=True,
    )
    return MPtoASEConverter(atoms=atoms).convert_input_generator(input_generator)


@job
def odac_static_job(
    atoms: Atoms,
    kpts: tuple = (1, 1, 1),
    copy_files: CopyFiles | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs: Any,
) -> VaspSchema:
    """
    Carry out a static calculation with ODAC settings.

    Parameters
    ----------
    atoms
        Atoms object
    kpts
        The k-point grid mesh. Please choose this carefully. The original
        ODAC23 paper uses a 1x1x1 k-point grid, which is not generally suitable.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """

    calc_defaults = {
        "kpts": kpts,
        "nwrite": 2,
        "xc": "pbe",
        "ivdw": 12,
        "encut": 600.0,
        "lcharg": False,
        "lwave": False,
        "ismear": 0,
        "sigma": 0.2,
        "ispin": 2,
        "prec": "accurate",
        "nelm": 120,
        "nelmin": 2,
        "ediff": 1e-5,
        "ediffg": -0.05,
        "maxmix": 40,
        "nsw": 0,
        "ibrion": 2,
        "isif": 3,
        "potim": 0.01,
        "algo": "normal",
        "ldiag": True,
        "lreal": "auto",
        "lplane": True,
        "gamma": True,
        "isym": 0,
        "pp_version": "54",
    }
    calc_defaults |= {"incar_copilot_mode": "critical", "use_custodian": False}
    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "ODAC Static"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
@requires(
    has_fairchem_oc,
    "fairchem-data-oc is not installed. Run `pip install quacc[fairchem]`",
)
def oc20_static_job(
    atoms: Atoms,
    copy_files: CopyFiles | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs: Any,
) -> VaspSchema:
    """
    Carry out a static calculation with OC20 settings.

    Parameters
    ----------
    atoms
        Atoms object
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """
    from fairchem.data.oc.utils.vasp import calculate_surface_k_points
    from fairchem.data.oc.utils.vasp_flags import VASP_FLAGS

    calc_defaults = VASP_FLAGS | {
        "kpts": calculate_surface_k_points(atoms),
        "xc": "RPBE",
        "pp_version": "54",
        "incar_copilot_mode": "critical",
        "use_custodian": False,
    }

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "OC20 Static"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
def oc22_static_job(
    atoms: Atoms,
    copy_files: CopyFiles | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs: Any,
) -> VaspSchema:
    """
    Carry out a static calculation with OC22 settings.

    Follows the [OC22 input generator](https://github.com/Open-Catalyst-Project/Open-Catalyst-Dataset/blob/23cdca3e7f14f6f5d8c3971ef2dbdcdba9dfeee3/ocdata/oc22_dataset/sets.py).
    Dipole corrections are enabled; pass ``ldipol=False`` for clean symmetric slabs.

    Parameters
    ----------
    atoms
        Oxide slab with its surface normal along z. Uses quacc's magnetic-moment
        handling with pymatgen element defaults.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """
    from pymatgen.io.vasp.sets import MVLSlabSet

    a, b, _ = atoms.cell.lengths()
    calc_defaults = MPtoASEConverter(atoms=atoms).convert_input_set(
        MVLSlabSet(
            set_mix=False,
            auto_dipole=True,
            user_potcar_functional="PBE_54",
            user_potcar_settings={"W": "W_sv"},
            user_incar_settings={
                "GGA": "PE",
                "ENCUT": 500,
                "EDIFF": 1e-4,
                "ISIF": 0,
                "NSW": 0,
                "SYMPREC": 1e-10,
                "NCORE": 4,  # Match the OC22 generator; MVLSlabSet has no NCORE default.
                "NELM": 60,
                "LREAL": False,
                "LASPH": False,
                "MAGMOM": None,
            },
        )
    )
    calc_defaults |= {
        "elemental_magmoms": MVLSlabSet.CONFIG["INCAR"]["MAGMOM"],
        "preset_mag_default": 0.6,
        "kpts": (ceil(30 / a), ceil(30 / b), 1),
        "incar_copilot_mode": "off",
        "use_custodian": False,
    }

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "OC22 Static"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
@requires(
    has_fairchem_oc,
    "fairchem-data-oc is not installed. Run `pip install quacc[fairchem]`",
)
def oc25_static_job(
    atoms: Atoms,
    copy_files: CopyFiles | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs: Any,
) -> VaspSchema:
    """
    Carry out a static calculation with OC25 evaluation settings.

    Uses ``ediff=1e-6`` for validation/test single points, following
    [Section 2.2.4 of the OC25 paper](https://arxiv.org/html/2509.17862v1#S2.SS2.SSS4).
    Pass ``ediff=1e-4`` to use the training-data electronic convergence threshold.

    Parameters
    ----------
    atoms
        Atoms object for the interface, with the surface normal along z.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.

    References
    ----------
    - [OC25 paper](https://arxiv.org/abs/2509.17862)
    - [OC25 input generation](https://github.com/facebookresearch/fairchem/issues/1616#issuecomment-3524048489)
    """
    from fairchem.data.oc.utils.vasp import calculate_surface_k_points
    from fairchem.data.oc.utils.vasp_flags import SOLVENT_BASE_FLAGS

    calc_defaults = SOLVENT_BASE_FLAGS | {
        "ediff": 1e-6,
        "nsw": 0,
        "kpts": calculate_surface_k_points(atoms),
        "xc": "RPBE",
        "pp_version": "64",
        "setups": "recommended",
        "incar_copilot_mode": "off",
        "use_custodian": False,
    }

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "OC25 Static"} | (additional_fields or {}),
        copy_files=copy_files,
    )
