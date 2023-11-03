"""Core common workflows"""
<<<<<<< Updated upstream
from typing import TYPE_CHECKING

from quacc import job

=======
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import job
from quacc.runners.ase import run_calc

>>>>>>> Stashed changes

if TYPE_CHECKING:
    from ase import Atoms
    from ase.calculators.calculator import Calculator
    from numpy.typing import NDArray

<<<<<<< Updated upstream
    from quacc.runners.ase import run_calc

=======
>>>>>>> Stashed changes

@job
def force_job(atoms: Atoms, calculator: Calculator) -> NDArray:
    """
    Calculate the forces.

    Parameters
    ----------
    atoms
        Atoms object
    calculator
        Calculator to use.

    Returns
    -------
    NDArray
        Forces
    """
    atoms.calc = calculator
    return run_calc(atoms).get_forces()
