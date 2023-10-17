"""Type hints for quacc.schemas.cclib"""
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import TypedDict

    from numpy.typing import Any, NDArray
    from pymatgen.core import Molecule

    from quacc.schemas._aliases.ase import parameters
    from quacc.schemas._aliases.atoms import AtomsSchema

    class AdditionalAttributes(TypedDict, total=False):
        energy: float
        homo_energies: list[float]
        lumo_energies: list[float]
        homo_lumo_gaps: list[float]
        min_homo_lumo_gap: float

    class CalcAttributes(TypedDict, total=False):
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

    class AllAttributes(Attributes, AdditionalAttributes, CalcAttributes):
        pass

    class cclibSchema(TypedDict, total=False):
        dir_name: str
        nid: str
        logfile: str
        molecule: Molecule
        molecule_initial: Molecule
        molecule_unoriented: Molecule
        parameters: parameters
        results: AllAttributes
        trajectory: list[AtomsSchema]
        task_label: str
        tags: list[str]
