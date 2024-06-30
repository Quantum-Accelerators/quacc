"""Custodian handlers for QChem."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import get_settings

if TYPE_CHECKING:
    from pathlib import Path


class _DefaultSettingType:
    pass


_DEFAULT_SETTING = _DefaultSettingType()
has_ob = bool(find_spec("openbabel"))


@requires(
    has_ob, "Openbabel must be installed. Try conda install -c conda-forge openbabel"
)
def run_custodian(
    qchem_cmd: str | _DefaultSettingType = _DEFAULT_SETTING,
    qchem_cores: int | _DefaultSettingType = _DEFAULT_SETTING,
    qchem_local_scratch: str | Path | _DefaultSettingType = _DEFAULT_SETTING,
    qchem_use_error_handlers: bool | _DefaultSettingType = _DEFAULT_SETTING,
    qchem_custodian_max_errors: int | _DefaultSettingType = _DEFAULT_SETTING,
    qchem_nbo_exe: str | Path | _DefaultSettingType = _DEFAULT_SETTING,
    directory: str | Path | None = None,
) -> list[list[dict]]:
    """
    Function to run QChem Custodian.

    Parameters
    ----------
    qchem_cmd
        Q-Chem command. Defaults to "qchem" in settings.
    qchem_cores
        Number of cores to use for the Q-Chem calculation.
    qchem_local_scratch
        Compute-node local scratch directory in which Q-Chem should perform IO.
        Defaults to /tmp in settings.
    qchem_use_error_handlers
        Whether or not to employ error handlers. Defaults to True in settings.
    qchem_custodian_max_errors
        Maximum number of errors to allow before stopping the run. Defaults to 5
        in settings.
    qchem_nbo_exe
        The full path to the NBO executable.
    directory
        The runtime directory.

    Returns
    -------
    list[list[dict]]
        Lists of lists of errors.
    """
    # Adapted from atomate.qchem.firetasks.run_calc
    from custodian import Custodian
    from custodian.qchem.handlers import QChemErrorHandler
    from custodian.qchem.jobs import QCJob

    settings = get_settings()

    # Set defaults
    qchem_cores = (
        settings.QCHEM_NUM_CORES if qchem_cores == _DEFAULT_SETTING else qchem_cores
    )
    qchem_cmd = settings.QCHEM_CMD if qchem_cmd == _DEFAULT_SETTING else qchem_cmd
    qchem_local_scratch = (
        settings.QCHEM_LOCAL_SCRATCH
        if qchem_local_scratch == _DEFAULT_SETTING
        else qchem_local_scratch
    )
    qchem_use_error_handlers = (
        settings.QCHEM_USE_ERROR_HANDLERS
        if qchem_use_error_handlers == _DEFAULT_SETTING
        else qchem_use_error_handlers
    )
    qchem_custodian_max_errors = (
        settings.QCHEM_CUSTODIAN_MAX_ERRORS
        if qchem_custodian_max_errors == _DEFAULT_SETTING
        else qchem_custodian_max_errors
    )
    qchem_nbo_exe = (
        settings.QCHEM_NBO_EXE if qchem_nbo_exe == _DEFAULT_SETTING else qchem_nbo_exe
    )

    # Error handlers for Q-Chem
    handlers = [QChemErrorHandler()] if qchem_use_error_handlers else []

    # Run Q-Chem
    jobs = [
        QCJob(
            qchem_command=qchem_cmd,
            max_cores=qchem_cores,
            calc_loc=str(qchem_local_scratch),
            nboexe=str(qchem_nbo_exe),
        )
    ]

    c = Custodian(
        handlers,
        jobs,
        max_errors=qchem_custodian_max_errors,
        terminate_on_nonzero_returncode=False,
        directory=directory,
    )

    return c.run()


if __name__ == "__main__":
    run_custodian()
