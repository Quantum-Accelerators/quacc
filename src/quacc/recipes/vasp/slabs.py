"""Recipes for slabs"""
from __future__ import annotations

import covalent as ct
from ase import Atoms

from quacc.calculators.vasp import Vasp
from quacc.schemas.vasp import VaspSchema, summarize_run
from quacc.util.calc import run_calc
from quacc.util.slabs import make_adsorbate_structures, make_max_slabs_from_bulk


@ct.electron
def slab_static_job(
    atoms: Atoms | dict,
    preset: str | None = None,
    calc_swaps: dict | None = None,
) -> VaspSchema:
    """
    Function to carry out a single-point calculation on a slab.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    preset
        Preset to use.
    calc_swaps
        dictionary of custom kwargs for the calculator.

    Returns
    -------
    VaspSchema
        Dictionary of results from quacc.schemas.vasp.summarize_run
    """
    atoms = atoms if isinstance(atoms, Atoms) else atoms["atoms"]
    calc_swaps = calc_swaps or {}

    defaults = {
        "auto_dipole": True,
        "ismear": -5,
        "laechg": True,
        "lcharg": True,
        "lvhar": True,
        "lwave": True,
        "nedos": 5001,
        "nsw": 0,
    }
    flags = defaults | calc_swaps

    calc = Vasp(atoms, preset=preset, **flags)
    atoms.calc = calc
    atoms = run_calc(atoms)

    return summarize_run(atoms, additional_fields={"name": "VASP Slab Static"})


@ct.electron
def slab_relax_job(
    atoms: Atoms | dict,
    preset: str | None = None,
    calc_swaps: dict | None = None,
) -> VaspSchema:
    """
    Function to relax a slab.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    preset
        Preset to use.
    calc_swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    VaspSchema
        Dictionary of results from quacc.schemas.vasp.summarize_run
    """
    atoms = atoms if isinstance(atoms, Atoms) else atoms["atoms"]
    calc_swaps = calc_swaps or {}

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
    flags = defaults | calc_swaps

    calc = Vasp(atoms, preset=preset, **flags)
    atoms.calc = calc
    atoms = run_calc(atoms)

    return summarize_run(atoms, additional_fields={"name": "VASP Slab Relax"})


def bulk_to_slabs_flow(
    atoms: Atoms | dict,
    slabgen_kwargs: dict | None = None,
    slab_relax: ct.electron = slab_relax_job,
    slab_static: ct.electron | None = slab_static_job,
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
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    slabgen_kwargs
        Additional keyword arguments to pass to make_max_slabs_from_bulk()
    slab_relax
        Default to use for the relaxation of the slab structures.
    slab_static
        Default to use for the static calculation of the slab structures.
    slab_relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    slab_static_kwargs
        Additional keyword arguments to pass to the static calculation.

    Returns
    -------
    list[VaspSchema]
        List of dictionary results from quacc.schemas.vasp.summarize_run
    """
    atoms = atoms if isinstance(atoms, Atoms) else atoms["atoms"]
    slab_relax_kwargs = slab_relax_kwargs or {}
    slab_static_kwargs = slab_static_kwargs or {}
    slabgen_kwargs = slabgen_kwargs or {}

    @ct.electron
    @ct.lattice
    def _relax_distributed(slabs):
        return [slab_relax(slab, **slab_relax_kwargs) for slab in slabs]

    @ct.electron
    @ct.lattice
    def _relax_and_static_distributed(slabs):
        return [
            slab_static(
                slab_relax(slab, **slab_relax_kwargs),
                **slab_static_kwargs,
            )
            for slab in slabs
        ]

    slabs = ct.electron(make_max_slabs_from_bulk)(atoms, **slabgen_kwargs)

    if slab_static is None:
        return _relax_distributed(slabs)
    return _relax_and_static_distributed(slabs)


def slab_to_ads_flow(
    slab: Atoms,
    adsorbate: Atoms,
    make_ads_kwargs: dict | None = None,
    slab_relax: ct.electron = slab_relax_job,
    slab_static: ct.electron | None = slab_static_job,
    slab_relax_kwargs: dict | None = None,
    slab_static_kwargs: dict | None = None,
) -> list[VaspSchema]:
    """
    Workflow consisting of:
    1. Slab-adsorbate generation
    2. Slab-adsorbate relaxations
    3. Slab-adsorbate statics (optional)

    Parameters
    ----------
    slab
        Atoms object for the slab structure.
    adsorbate
        Atoms object for the adsorbate.
    make_ads_kwargs
        Additional keyword arguments to pass to make_adsorbate_structures()
    slab_relax
        Default to use for the relaxation of the slab structure.
    slab_static
        Default to use for the static calculation of the slab structures.
    slab_relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    slab_static_kwargs
        Additional keyword arguments to pass to the static calculation.

    Returns
    -------
    list[VaspSchema]
        List of dictionaries of results from quacc.schemas.vasp.summarize_run
    """

    slab_relax_kwargs = slab_relax_kwargs or {}
    slab_static_kwargs = slab_static_kwargs or {}
    make_ads_kwargs = make_ads_kwargs or {}

    @ct.electron
    @ct.lattice
    def _relax_distributed(slabs):
        return [slab_relax(slab, **slab_relax_kwargs) for slab in slabs]

    @ct.electron
    @ct.lattice
    def _relax_and_static_distributed(slabs):
        return [
            slab_static(
                slab_relax(slab, **slab_relax_kwargs),
                **slab_static_kwargs,
            )
            for slab in slabs
        ]

    ads_slabs = ct.electron(make_adsorbate_structures)(
        slab, adsorbate, **make_ads_kwargs
    )

    if slab_static is None:
        return _relax_distributed(ads_slabs)
    return _relax_and_static_distributed(ads_slabs)
