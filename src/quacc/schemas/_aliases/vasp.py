"""Type hints for quacc.schemas.vasp."""
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
