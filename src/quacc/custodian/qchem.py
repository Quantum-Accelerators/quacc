"""Custodian handlers for QChem"""
from __future__ import annotations

import sys

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
    qchem_cores: int = 1,
    qchem_cmd: str | None = None,
    qchem_local_scratch: str | None = None,
    qchem_use_error_handlers: bool | None = None,
    qchem_custodian_max_errors: int | None = None,
) -> None:
    """
    Function to run QChem Custodian

    Parameters
    ----------
    qchem_cores
        Number of cores to use for the Q-Chem calculation.
    qchem_cmd
        Q-Chem command. Defaults to "qchem" in settings.
    qchem_local_scratch
        Compute-node local scratch directory in which Q-Chem should perform IO. Defaults to /tmp in settings.
    qchem_use_error_handlers
        Whether or not to employ error handlers. Defaults to True in settings.
    qchem_custodian_max_errors
        Maximum number of errors to allow before stopping the run. Defaults to 5 in settings.

    Returns
    -------
    None
    """
    # Adapted from atomate.qchem.firetasks.run_calc
    from custodian import Custodian
    from custodian.qchem.handlers import QChemErrorHandler
    from custodian.qchem.jobs import QCJob

    # Set defaults
    qchem_cmd = SETTINGS.QCHEM_CMD if qchem_cmd is None else qchem_cmd
    qchem_local_scratch = (
        SETTINGS.QCHEM_LOCAL_SCRATCH
        if qchem_local_scratch is None
        else qchem_local_scratch
    )
    qchem_use_error_handlers = (
        SETTINGS.QCHEM_USE_ERROR_HANDLERS
        if qchem_use_error_handlers is None
        else qchem_use_error_handlers
    )
    qchem_custodian_max_errors = (
        SETTINGS.QCHEM_CUSTODIAN_MAX_ERRORS
        if qchem_custodian_max_errors is None
        else qchem_custodian_max_errors
    )

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
            calc_loc=qchem_local_scratch,
        )
    ]

    c = Custodian(
        handlers,
        jobs,
        max_errors=qchem_custodian_max_errors,
    )

    c.run()


if __name__ == "__main__":
    run_custodian(qchem_cores=sys.argv[1]) if len(sys.argv) > 1 else run_custodian()
