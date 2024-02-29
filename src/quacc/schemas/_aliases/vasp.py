"""Type hints for quacc.schemas.vasp."""

from __future__ import annotations

from typing import TypedDict

from quacc.schemas._aliases.ase import RunSchema
from quacc.schemas._aliases.emmet import TaskDoc


class BaderSchema(TypedDict, total=False):
    """Type hint associated with quacc.schemas.vasp._bader_runner."""

    atomic_volume: float
    partial_charges: list[float]
    spin_moments: list[float]
    bader_version: float
    min_dist: list[float]


class DDECSchema(TypedDict, total=False):
    """Type hint associated with quacc.schemas.vasp._ddec_runner."""

    partial_charges: list[float]
    spin_moments: list[float]
    dipoles: list[float]
    bond_order_sums: list[float]
    bond_order_dict: dict
    rsquared_moments: list[float]
    rcubed_moments: list[float]
    rfourth_moments: list[float]


class CM5Schema(TypedDict):
    """Type hint used in DDECSchema"""

    partial_charges: list[float]


class ChargemolSchema(TypedDict, total=False):
    """Type hint associated with quacc.schemas.vasp._chargemol_runner`"""

    ddec: DDECSchema
    cm5: CM5Schema


class VaspSchema(RunSchema, TaskDoc):
    """Type hint associated with [quacc.schemas.vasp.vasp_summarize_run][]"""

    bader: BaderSchema
    chargemol: ChargemolSchema


class MPGGARelaxFlowSchema(VaspSchema):
    """Type hint associated with the MP GGA relaxation flows."""

    relax1: VaspSchema
    relax2: VaspSchema
    static: VaspSchema


class MPMetaGGARelaxFlowSchema(MPGGARelaxFlowSchema):
    """Type hint associated with the MP meta-GGA relaxation flows."""

    prerelax: VaspSchema


class QMOFRelaxSchema(VaspSchema):
    """Type hint associated with the QMOF relaxation jobs."""

    prerelax_lowacc: VaspSchema | None
    position_relax_lowacc: VaspSchema
    volume_relax_lowacc: VaspSchema | None
    double_relax: VaspSchema
