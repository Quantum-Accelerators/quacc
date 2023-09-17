"""Recipes for slabs"""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import flow, job, subflow
from quacc.recipes.vasp.core import _base_job
from quacc.schemas import fetch_atoms
from quacc.utils.slabs import make_adsorbate_structures, make_slabs_from_bulk

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.vasp import VaspSchema


@job
def slab_static_job(
    atoms: Atoms | dict,
    preset: str | None = None,
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> VaspSchema:
    """
    Function to carry out a single-point calculation on a slab.

    ??? Note

        Calculator Defaults:

        ```python
        {
            "auto_dipole": True,
            "ismear": -5,
            "laechg": True,
            "lcharg": True,
            "lreal": False,
            "lvhar": True,
            "lwave": True,
            "nedos": 5001,
            "nsw": 0,
        }
        ```

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    preset
        Preset to use.
    calc_swaps
        dictionary of custom kwargs for the calculator.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][]
    """

    defaults = {
        "auto_dipole": True,
        "ismear": -5,
        "laechg": True,
        "lcharg": True,
        "lreal": False,
        "lvhar": True,
        "lwave": True,
        "nedos": 5001,
        "nsw": 0,
    }
    return _base_job(
        atoms,
        preset=preset,
        defaults=defaults,
        calc_swaps=calc_swaps,
        additional_fields={"name": "VASP Slab Static"},
        copy_files=copy_files,
    )


@job
def slab_relax_job(
    atoms: Atoms | dict,
    preset: str | None = None,
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> VaspSchema:
    """
    Function to relax a slab.

    ??? Note

        Calculator Parameters:

        ```python
        {
            "auto_dipole": True,
            "ediffg": -0.02,
            "isif": 2,
            "ibrion": 2,
            "isym": 0,
            "lcharg": False,
            "lwave": False,
            "nsw": 200,
        }
        ```

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    preset
        Preset to use.
    calc_swaps
        Dictionary of custom kwargs for the calculator.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][]
    """

    defaults = {
        "auto_dipole": True,
        "ediffg": -0.02,
        "isif": 2,
        "ibrion": 2,
        "isym": 0,
        "lcharg": False,
        "lwave": False,
        "nsw": 200,
    }
    return _base_job(
        atoms,
        preset=preset,
        defaults=defaults,
        calc_swaps=calc_swaps,
        additional_fields={"name": "VASP Slab Relax"},
        copy_files=copy_files,
    )


@flow
def bulk_to_slabs_flow(
    atoms: Atoms | dict,
    make_slabs_kwargs: dict | None = None,
    run_static: bool = True,
    slab_relax_kwargs: dict | None = None,
    slab_static_kwargs: dict | None = None,
) -> list[VaspSchema]:
    """
    Workflow consisting of:

    1. Slab generation

    2. Slab relaxations

    3. Slab statics (optional)

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    make_slabs_kwargs
        Additional keyword arguments to pass to [quacc.utils.slabs.make_slabs_from_bulk][]
    run_static
        Whether to run the static calculation.
    slab_relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    slab_static_kwargs
        Additional keyword arguments to pass to the static calculation.

    Returns
    -------
    list[VaspSchema]
        List of dictionary results from [quacc.schemas.vasp.vasp_summarize_run][]
    """
    slab_relax_kwargs = slab_relax_kwargs or {}
    slab_static_kwargs = slab_static_kwargs or {}
    make_slabs_kwargs = make_slabs_kwargs or {}

    @job
    def _make_slabs(atoms):
        atoms = fetch_atoms(atoms)
        return make_slabs_from_bulk(atoms, **make_slabs_kwargs)

    @subflow
    def _relax_distributed(slabs):
        return [slab_relax_job(slab, **slab_relax_kwargs) for slab in slabs]

    @subflow
    def _relax_and_static_distributed(slabs):
        return [
            slab_static_job(
                slab_relax_job(slab, **slab_relax_kwargs),
                **slab_static_kwargs,
            )
            for slab in slabs
        ]

    slabs = _make_slabs(atoms)

    return (
        _relax_and_static_distributed(slabs)
        if run_static
        else _relax_distributed(slabs)
    )


@flow
def slab_to_ads_flow(
    slab: Atoms,
    adsorbate: Atoms,
    make_ads_kwargs: dict | None = None,
    run_static: bool = True,
    slab_relax_kwargs: dict | None = None,
    slab_static_kwargs: dict | None = None,
) -> list[VaspSchema]:
    """
    Workflow consisting of: 1. Slab-adsorbate generation 2. Slab-adsorbate
    relaxations 3. Slab-adsorbate statics (optional)

    Parameters
    ----------
    slab
        Atoms object for the slab structure.
    adsorbate
        Atoms object for the adsorbate.
    make_ads_kwargs
        Additional keyword arguments to pass to [quacc.utils.slabs.make_adsorbate_structures][]
    run_static
        Whether to run the static calculation.
    slab_relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    slab_static_kwargs
        Additional keyword arguments to pass to the static calculation.

    Returns
    -------
    list[VaspSchema]
        List of dictionaries of results from [quacc.schemas.vasp.vasp_summarize_run][]
    """

    slab_relax_kwargs = slab_relax_kwargs or {}
    slab_static_kwargs = slab_static_kwargs or {}
    make_ads_kwargs = make_ads_kwargs or {}

    @job
    def _make_ads_slabs(atoms, adsorbate):
        atoms = fetch_atoms(atoms)
        return make_adsorbate_structures(atoms, adsorbate, **make_ads_kwargs)

    @subflow
    def _relax_distributed(slabs):
        return [slab_relax_job(slab, **slab_relax_kwargs) for slab in slabs]

    @subflow
    def _relax_and_static_distributed(slabs):
        return [
            slab_static_job(
                slab_relax_job(slab, **slab_relax_kwargs),
                **slab_static_kwargs,
            )
            for slab in slabs
        ]

    ads_slabs = _make_ads_slabs(slab, adsorbate)

    return (
        _relax_and_static_distributed(ads_slabs)
        if run_static
        else _relax_distributed(ads_slabs)
    )
