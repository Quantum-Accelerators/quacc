"""Recipes to reproduce OMol settings"""

from fairchem.data.omol.orca.recipes import single_point_calculation
from typing import TYPE_CHECKING
from importlib.util import find_spec

if TYPE_CHECKING:
  from ase.atoms import Atoms
  from fairchem.data.omol.orca.calc import Vertical
  
has_fairchem = bool(find_spec("fairchem"))
has_fairchem_omol = (
    has_fairchem
    and bool(find_spec("fairchem.data"))
    and bool(find_spec("fairchem.data.omol"))
)

@job
@requires(
    has_fairchem_omol,
    "fairchem-data-omol is not installed. Run `pip install quacc[fairchem]`",
)
def omol_static_job(
  atoms: Atoms,
  charge: int = 0,
  spin_multiplicity: int = 1,
  vertical: Vertical = Vertical.Default,
  orcasimpleinput: list[str] | None = None,
  orcablocks: list[str] | None = None,
  copy_files: SourceDirectory | Copy | None = None,
  additional_fields: dict[str, Any] | None = None,
)

  return single_point_calculation(atoms, charge, spin_multiplicity, vertical=vertical, orcasimpleinput=orcasimpleinput, orcablocks=orcablocks,copy_files=copy_files,additional_fields=additional_fields)
