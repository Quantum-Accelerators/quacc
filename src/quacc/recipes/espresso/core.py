"""Core recipes for espresso."""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase import Atoms

from quacc import SETTINGS, job
from quacc.calculators.espresso.espresso import (Espresso, EspressoProfile,
                                                 EspressoTemplate)
from quacc.runners.ase import run_calc
from quacc.schemas.ase import summarize_run

if TYPE_CHECKING:
    from typing import Any

    from quacc.schemas._aliases.ase import RunSchema


@job
def static_job(
    atoms: Atoms,
    preset: str | None = None,
    copy_files: list[str] | None = None,
    tmpdir: str | None = None,
    calc_clean: bool = True,
    **kwargs,
) -> RunSchema:
    """
    Function to carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    **kwargs
        Custom kwargs for the espresso calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.espresso.Espresso` calculator.

        !!! Info "Calculator defaults"

            ```python
            {
                "ecutwfc": 40,
                "ecutrho": 160,
            }
            ```

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    calc_defaults = {
        'input_data': {
            'control': {
                'calculation': 'scf',
                'restart_mode': 'from_scratch',
                'pseudo_dir': SETTINGS.ESPRESSO_PP_PATH,
            },
            'system': {
                "ecutwfc": 60,
                "ecutrho": 240,
            },
            'electrons': {
                'conv_thr': 1e-8,
                'mixing_mode': 'plain',
                'mixing_beta': 0.7
            }
        }
    }
    
    template = EspressoTemplate('pw')

    return _base_job(
        atoms,
        preset=preset,
        template=template,
        calc_defaults=calc_defaults,
        calc_swaps=kwargs,
        additional_fields={"name": "pw.x static"},
        copy_files=copy_files,
        tmpdir=tmpdir,
        calc_clean=calc_clean
    )

@job
def ph_job(
    preset: str | None = None,
    copy_files: list[str] | None = None,
    tmpdir: str | None = None,
    calc_clean: bool = True,
    **kwargs,
) -> RunSchema:
    """
    Function to carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    **kwargs
        Custom kwargs for the espresso calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.espresso.Espresso` calculator.

        !!! Info "Calculator defaults"

            ```python
            {
                "ecutwfc": 40,
                "ecutrho": 160,
            }
            ```

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    # Default is phonon at gamma, with tight convergence
    # and a lower alphamix(1) as it is very very often
    # recommended in the QE mailing list...
    calc_defaults = {
        'input_data': {
            'inputph': {
                'tr2_ph': 1e-16,
                'alpha_mix(1)': 0.1,
                'nmix_ph': 12,
                'verbosity': 'high',

            }
        },
        'qpts': [(0, 0, 0)]
    }
    
    template = EspressoTemplate('ph')
    profile = EspressoProfile(
        binary=str(SETTINGS.ESPRESSO_PH_CMD),
        pseudo_path=None
        )

    return _base_job(
        Atoms(),
        preset = preset,
        template = template,
        profile = profile,
        calc_defaults=calc_defaults,
        calc_swaps=kwargs,
        additional_fields={"name": "ph.x static"},
        copy_files=copy_files,
        tmpdir=tmpdir,
        calc_clean = calc_clean
    )

def _base_job(
    atoms: Atoms,
    preset: str | None = None,
    template: EspressoTemplate | None = None,
    profile: EspressoProfile | None = None,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
    tmpdir: str | None = None,
    calc_clean: bool = True
) -> RunSchema:
    """
    Base function to carry out espresso recipes.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Name of the preset to use
    calc_defaults
        The default calculator parameters.
    calc_swaps
        Custom kwargs for the espresso calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.espresso.Espresso` calculator.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    atoms.calc = Espresso(input_atoms = atoms,
                          preset = preset,
                          template = template,
                          profile = profile,
                          calc_defaults = calc_defaults,
                          **calc_swaps)
    
    final_atoms = run_calc(atoms,
                           copy_files=copy_files,
                           tmpdir=tmpdir,
                           calc_clean=calc_clean)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields=additional_fields,
    )
