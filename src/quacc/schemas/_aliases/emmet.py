"""Aliases for type hinting for emmet."""

from __future__ import annotations

from datetime import datetime
from typing import Any, Literal, TypedDict

from emmet.core.math import ListMatrix3D, Matrix3D, Vector3D
from emmet.core.symmetry import CrystalSystem
from emmet.core.vasp.calc_types import CalcType
from emmet.core.vasp.calc_types.enums import RunType, TaskType
from emmet.core.vasp.calculation import VaspObject
from emmet.core.vasp.task_valid import TaskState
from pymatgen.core.composition import Composition
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, Potcar


class SymmetryData(TypedDict):
    """Type hint associated with [emmet.core.symmetry.SymmetryData][]"""

    crystal_system: CrystalSystem
    symbol: str
    number: int
    point_group: str
    symprec: float
    version: str


class PointGroupData(TypedDict):
    """Type hint associated with [emmet.core.symmetry.PointGroupData][]"""

    point_group: str
    rotation_number: float
    linear: bool
    tolerance: float
    eigen_tolerance: float
    matrix_tolerance: float


class EmmetBase(TypedDict):
    """Type hint associated with `emmet.core.base.EmmetBaseModel`."""

    emmet_version: str
    pymatgen_version: str
    pull_request: int | None
    database_version: str | None
    build_date: datetime
    license: Literal["BY-C", "BY-NC"]


class StructureMetadata(EmmetBase):
    """Type hint associated with [emmet.core.structure.StructureMetadata][]"""

    nsites: int
    elements: list[Element]
    nelements: int
    composition: Composition
    formula_pretty: str
    formula_anonymous: str
    chemsys: str
    volume: float
    density: float
    density_atomic: float
    symmetry: SymmetryData


class MoleculeMetadata(EmmetBase):
    """Type hint associated with [emmet.core.structure.MoleculeMetadata][]"""

    charge: int
    spin_multiplicity: int
    natoms: int
    elements: list[Element]
    nelements: int
    nelectrons: int
    composition: Composition
    composition_reduced: Composition
    formula_alphabetical: str
    formula_pretty: str
    formula_anonymous: str
    chemsys: str
    symmetry: PointGroupData


class PotcarSpec(TypedDict):
    """Type hint associated with emmet.core.vasp.calculation.PotcarSpec."""

    titel: str
    hash: str


class CalculationInput(TypedDict):
    """Type hint associated with emmet.core.vasp.calculation.CalculationInput."""

    incar: dict[str, Any]
    kpoints: dict[str, Any]
    nkpoints: int
    potcar: list[str]
    potcar_spec: list[PotcarSpec]
    potcar_type: list[str]
    parameters: dict[str, Any]
    lattice_rec: Lattice
    structure: Structure
    is_hubbard: bool
    hubbards: dict[str, float]


class ElectronicStep(TypedDict):
    """Type hint associated with emmet.core.vasp.calculation.ElectronicStep."""

    alphaZ: float
    ewald: float
    hartreedc: float
    XCdc: float
    pawpsdc: float
    pawaedc: float
    eentropy: float
    bandstr: float
    atom: float
    e_fr_energy: float
    e_wo_entrp: float
    e_0_energy: float


class IonicStep(TypedDict):
    """Type hint associated with emmet.core.vasp.calculation.IonicStep."""

    e_fr_energy: float
    e_wo_entrp: float
    e_0_energy: float
    forces: list[Vector3D]
    stress: Matrix3D
    electronic_steps: list[ElectronicStep]
    structure: Structure


class FrequencyDependentDielectric(TypedDict):
    """Type hint associated with
    emmet.core.vasp.calculation.FrequencyDependentDielectric."""

    real: list[list[float]]
    imaginary: list[list[float]]
    energy: list[float]


class ElectronPhononDisplacedStructures(TypedDict):
    """Type hint associated with
    emmet.core.vasp.calculation.ElectronPhononDisplacedStructures."""

    temperatures: list[float]
    structures: list[Structure]


class RunStatistics(TypedDict):
    """Type hint associated with emmet.core.vasp.calculation.RunStatistics."""

    average_memory: float
    max_memory: float
    elapsed_time: float
    system_time: float
    user_time: float
    total_time: float
    cores: int


class CalculationOutput(TypedDict, total=False):
    """Type hint associated with emmet.core.vasp.calculation.CalculationOutput."""

    energy: float
    energy_per_atom: float
    structure: Structure
    efermi: float
    is_metal: bool
    bandgap: float
    cbm: float
    vbm: float
    is_gap_direct: bool
    direct_gap: float
    transition: str
    mag_density: float
    epsilon_static: ListMatrix3D
    epsilon_static_wolfe: ListMatrix3D
    epsilon_ionic: ListMatrix3D
    frequency_dependent_dielectric: FrequencyDependentDielectric
    ionic_steps: list[IonicStep]
    locpot: dict[int, list[float]]
    outcar: dict[str, Any]
    force_constants: list[list[Matrix3D]]
    normalmode_frequencies: list[float]
    normalmode_eigenvals: list[float]
    normalmode_eigenvecs: list[list[Vector3D]]
    elph_displaced_structures: ElectronPhononDisplacedStructures
    dos_properties: dict[str, dict[str, dict[str, float]]]
    run_stats: RunStatistics


class Calculation(TypedDict):
    """Type hint associated with emmet.core.vasp.calculation.Calculation."""

    dir_name: str
    vasp_version: str
    has_vasp_completed: TaskState | bool
    input: CalculationInput
    output: CalculationOutput
    completed_at: str
    task_name: str
    output_file_paths: dict[str, str]
    bader: dict
    run_type: RunType
    task_type: TaskType
    calc_type: CalcType


class OrigInputs(TypedDict):
    """Type hint associated with emmet.core.tasks.OrigInputs."""

    incar: Incar
    poscar: Poscar
    kpoints: Kpoints
    potcar: Potcar


class OutputDoc(TypedDict):
    """Type hint associated with emmet.core.tasks.OutputDoc."""

    structure: Structure
    density: float
    energy: float
    forces: list[list[float]]
    stress: list[list[float]]
    energy_per_atom: float
    bandgap: float


class CustodianDoc(TypedDict):
    """Type hint associated with emmet.core.tasks.CustodianDoc."""

    corrections: list[Any]
    job: dict


class AnalysisDoc(TypedDict):
    """Type hint associated with emmet.core.tasks.AnalysisDoc."""

    delta_volume: float
    delta_volume_percent: float
    max_force: float
    warnings: list[str]
    errors: list[str]


class TaskDoc(StructureMetadata):
    """Type hint associated with emmet.core.tasks.TaskDoc."""

    tags: list[str] | None
    dir_name: str
    state: TaskState
    calcs_reversed: list[Calculation]
    structure: Structure
    task_type: CalcType | TaskType
    task_id: str
    orig_inputs: OrigInputs
    output: OutputDoc
    included_objects: VaspObject
    vasp_objects: dict[VaspObject, Any]
    entry: ComputedEntry
    task_label: str
    author: str
    transformations: dict[str, Any]
    additional_json: dict[str, Any]
    custodian: list[CustodianDoc]
    analysis: AnalysisDoc
    last_updated: datetime
