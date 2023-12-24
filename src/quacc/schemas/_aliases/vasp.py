"""Type hints for quacc.schemas.vasp."""
from __future__ import annotations

from typing import TypedDict

from quacc.schemas._aliases.ase import RunSchema
from quacc.schemas._aliases.emmet import TaskDoc


class BaderSchema(TypedDict, total=False):
    """Type hint associated with quacc.schemas.vasp._bader_runner."""

    atomic_volume: float
    bader_charge: float
    bader_spin: float
    bader_version: float
    min_dist: list[float]
    partial_charges: list[float]
    spin_moments: list[float]


class DDECSchema(TypedDict, total=False):
    partial_charges: list[float]
    spin_moments: list[float]
    dipoles: list[float]
    bond_order_sums: list[float]
    bond_order_dict: dict
    rsquared_moments: list[float]
    rcubed_moments: list[float]
    rfourth_moments: list[float]


class CM5Schema(TypedDict):
    partial_charges: list[float]


class ChargemolSchema(DDECSchema, CM5Schema):
    """Type hint associated with quacc.schemas.chargemol._chargemol_runner."""


class VaspSchema(RunSchema, TaskDoc):
    """Type hint associated with quacc.schemas.vasp.vasp_summarize_run."""

    bader: BaderSchema
    chargemol: ChargemolSchema


class DoubleRelaxSchema(VaspSchema):
    """Type hint associated with double relaxation jobs."""

    relax1: VaspSchema


class MPRelaxFlowSchema(VaspSchema):
    """Type hint associated with the MP relaxation flows"""

    prerelax: VaspSchema


class QMOFRelaxSchema(VaspSchema):
    """Type hint associated with the QMOF relaxation jobs"""

    prerelax_lowacc: VaspSchema | None
    position_relax_lowacc: VaspSchema
    volume_relax_lowacc: VaspSchema | None
    double_relax: VaspSchema
