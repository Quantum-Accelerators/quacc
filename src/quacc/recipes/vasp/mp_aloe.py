"""MP-ALOE-compatible VASP recipes."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import job
from quacc.calculators.vasp.params import MPtoASEConverter
from quacc.recipes.vasp._base import run_and_summarize

has_atomate2 = bool(find_spec("atomate2"))

if TYPE_CHECKING:
    from ase.atoms import Atoms

    from quacc.types import SourceDirectory, VaspSchema


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[mp]`")
def mp_aloe_static_job(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs
) -> VaspSchema:
    """
    Run a static calculation with MP-ALOE settings.

    The recipe uses the pymatgen MP24 relaxation settings with a plane-wave
    cutoff of 680 eV and ``KSPACING = 0.2``.

    Parameters
    ----------
    atoms
        Atoms object.
    prev_dir
        A previous directory for a prior step in the workflow.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to ``None`` to
        remove a pre-existing key entirely. User values take precedence over
        the MP-ALOE defaults.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
    """
    from atomate2.vasp.jobs.mp import MP24RelaxMaker

    calc_defaults = MPtoASEConverter(atoms=atoms, prev_dir=prev_dir).convert_maker(
        MP24RelaxMaker()
    )
    calc_defaults |= {
        "encut": 680,
        "kspacing": 0.2,
        "nsw": 0,
        "incar_copilot": "critical",
    }

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        report_mp_corrections=True,
        additional_fields={"name": "MP-ALOE Static"},
        copy_files={prev_dir: ["CHGCAR*", "WAVECAR*"]} if prev_dir else None,
    )
