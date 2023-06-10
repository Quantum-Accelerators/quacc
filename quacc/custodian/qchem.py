"""
Custodian handlers for QChem
"""
from __future__ import annotations

from custodian import Custodian
from monty.dev import requires

from quacc import SETTINGS

try:
    import openbabel as ob
except:
    ob = None


@requires(
    ob,
    "Openbabel must be installed. Try conda install -c conda-forge openbabel",
)
def run_custodian(
    qchem_cmd: str = SETTINGS.QChem_CMD,
    qchem_max_cores: int = SETTINGS.QChem_MAX_CORES,
    qchem_calc_loc: str = SETTINGS.QChem_CALC_LOC,
    qchem_custodian_max_errors: int = SETTINGS.QChem_CUSTODIAN_MAX_ERRORS,
    qchem_custodian_handlers: list[str] = SETTINGS.QChem_CUSTODIAN_HANDLERS,
) -> None:
    """
    Function to run QChem Custodian

    Parameters
    ----------
    qchem_cmd
        Q-Chem command. Defaults to "qchem" in settings.
    qchem_max_cores
        Maximum number of cores to use for the Q-Chem calculation. Defaults to 32 in settings.
    qchem_calc_loc
        Compute-node local scratch directory in which Q-Chem should perform IO. Defaults to /tmp in settings.
    qchem_custodian_max_errors
        Maximum number of errors to allow before stopping the run. Defaults to 5 in settings.
    qchem_custodian_handlers
        List of handlers to use in Custodian. See settings for list.

    Returns
    -------
    None
    """
    # Adapted from atomate.qchem.firetasks.run_calc

    from custodian.qchem.handlers import QChemErrorHandler
    from custodian.qchem.jobs import QCJob

    # Handlers for Q-Chem
    handlers_dict = {
        "QChemErrorHandler": QChemErrorHandler(),
    }

    handlers = []
    for handler_flag in qchem_custodian_handlers:
        if handler_flag not in handlers_dict:
            raise ValueError(f"Unknown Q-Chem error handler: {handler_flag}")
        handlers.append(handlers_dict[handler_flag])

    # Run Q-Chem

    jobs = [
        QCJob(
            qchem_command=qchem_cmd,
            max_cores=qchem_max_cores,
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
    run_custodian()
