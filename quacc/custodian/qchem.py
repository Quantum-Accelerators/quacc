"""
Custodian handlers for QChem
"""
from __future__ import annotations
import sys
import multiprocessing

from custodian import Custodian
from monty.dev import requires

from quacc import SETTINGS

try:
    import openbabel as ob
except ImportError:
    ob = None


@requires(
    ob,
    "Openbabel must be installed. Try conda install -c conda-forge openbabel",
)
def run_custodian(
    qchem_cores: int = multiprocessing.cpu_count(),
    qchem_cmd: str = SETTINGS.QCHEM_CMD,
    qchem_calc_loc: str = SETTINGS.QCHEM_CALC_LOC,
    qchem_use_error_handlers: bool = SETTINGS.QCHEM_USE_ERROR_HANDLERS,
    qchem_custodian_max_errors: int = SETTINGS.QCHEM_CUSTODIAN_MAX_ERRORS,
) -> None:
    """
    Function to run QChem Custodian

    Parameters
    ----------
    qchem_cores
        Number of cores to use for the Q-Chem calculation. Must be set by user on the command line.
    qchem_cmd
        Q-Chem command. Defaults to "qchem" in settings.
    qchem_calc_loc
        Compute-node local scratch directory in which Q-Chem should perform IO. Defaults to /tmp in settings.
    qchem_custodian_max_errors
        Maximum number of errors to allow before stopping the run. Defaults to 5 in settings.

    Returns
    -------
    None
    """
    # Adapted from atomate.qchem.firetasks.run_calc

    from custodian.qchem.handlers import QChemErrorHandler
    from custodian.qchem.jobs import QCJob

    # Error handlers for Q-Chem
    if qchem_use_error_handlers:
        handlers = [QChemErrorHandler()]
    else:
        handlers = []

    # Run Q-Chem
    jobs = [
        QCJob(
            qchem_command=qchem_cmd,
            max_cores=qchem_cores,
            calc_loc=qchem_calc_loc,
        )
    ]

    c = Custodian(
        handlers,
        jobs,
        max_errors=qchem_custodian_max_errors,
    )

    c.run()


if __name__ == "__main__":
    qchem_cores = sys.argv[1]
    run_custodian(qchem_cores=qchem_cores)
