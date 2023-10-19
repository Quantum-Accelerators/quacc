"""Type hints for quacc.schemas.cclib"""
from typing import Any, TypedDict

from numpy.typing import NDArray

from quacc.schemas._aliases.ase import RunSchema
from quacc.schemas._aliases.atoms import AtomsSchema


class AdditionalAttributes(TypedDict, total=False):
    """
    Additional type hints we custom-made based on cclib attributes.

    Uses cclib units.
    """

    final_scf_energy: float
    homo_energies: list[float] | None
    lumo_energies: list[float] | None
    homo_lumo_gaps: list[float] | None
    min_homo_lumo_gap: float | None


class PopAnalysisAttributes(TypedDict, total=False):
    """
    Type hints associated with cclib population analysis attribubtes.
    """

    aoresults: Any
    fragresults: Any
    fragcharges: Any
    density: Any
    donations: Any
    bdonations: Any
    repulsions: Any
    matches: Any
    refcharges: Any


class Attributes(TypedDict, total=False):
    """
    Type hints associated with cclib attribubtes.

    Refer to https://cclib.github.io/data.html
    """

    aonames: list[str]
    aooverlaps: NDArray
    atombasis: list[list[int]]
    atomcharges: dict[NDArray]
    atomcoords: NDArray
    atommasses: NDArray
    atomnos: NDArray
    atomspins: dict[NDArray]
    ccenergies: NDArray
    charge: int
    coreelectrons: NDArray
    dispersionenergies: NDArray
    enthalpy: float
    entropy: float
    etenergies: NDArray
    etoscs: NDArray
    etdips: NDArray
    etveldips: NDArray
    etmagdips: NDArray
    etrotats: NDArray
    etsecs: list[list]
    etsyms: list[str]
    freenergy: float
    fonames: list[str]
    fooverlaps: NDArray
    fragnames: list[str]
    frags: list[list[int]]
    gbasis: Any
    geotargets: NDArray
    geovalues: NDArray
    grads: NDArray
    hessian: NDArray
    homos: NDArray
    metadata: dict
    mocoeffs: list[NDArray]
    moenergies: list[NDArray]
    moments: list[NDArray]
    mosyms: list[list]
    mpenergies: NDArray
    mult: int
    natom: int
    nbasis: int
    nmo: int
    nmrtensors: dict[dict[NDArray]]
    nmrcouplingtensors: dict[dict[NDArray]]
    nocoeffs: NDArray
    nooccnos: NDArray
    nsocoeffs: list[NDArray]
    nsooccnos: list[NDArray]
    optdone: bool
    optstatus: NDArray
    polarizabilities: list[NDArray]
    pressure: float
    rotconsts: NDArray
    scancoords: NDArray
    scanenergies: list[float]
    scannames: list[str]
    scanparm: list[tuple]
    scfenergies: NDArray
    scftargets: NDArray
    scfvalues: list[NDArray]
    temperature: float
    time: NDArray
    transprop: Any
    vibanharms: NDArray
    vibdisps: NDArray
    vibfreqs: NDArray
    vibfconsts: NDArray
    vibirs: NDArray
    vibramans: NDArray
    vibrmasses: NDArray
    vibsyms: list[str]
    zpve: float


class AllAttributes(Attributes, AdditionalAttributes):
    """
    Type hint of all cclib attributes
    """


class cclibBaseSchema(TypedDict):
    """
    Type hint associated with quacc.schemas.cclib._make_cclib_schema
    """

    logfile: str
    attributes: AllAttributes
    pop_analysis: PopAnalysisAttributes | None
    trajectory: list[AtomsSchema]


class cclibSchema(cclibBaseSchema, RunSchema):
    """
    Type hint associated with quacc.schemas.cclib.cclib_summarize_run
    """
