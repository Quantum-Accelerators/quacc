"""Core recipes for VASP."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc.calculators.vasp import Vasp
from quacc.runners.ase import run_calc, run_opt
from quacc.schemas.vasp import summarize_vasp_opt_run, vasp_summarize_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.vasp import VaspASESchema, VaspSchema
    from quacc.utils.files import Filenames, SourceDirectory


class RunAndSummarize:
    """Run and summarize a VASP run."""

    def __init__(
        self,
        input_atoms: Atoms,
        preset: str | None = None,
        calc_defaults: dict[str, Any] | None = None,
        calc_swaps: dict[str, Any] | None = None,
        report_mp_corrections: bool = False,
        additional_fields: dict[str, Any] | None = None,
        copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    ) -> None:
        """
        Initialize the class.

        Parameters
        ----------
        input_atoms
            Atoms object
        preset
            Preset to use from `quacc.calculators.vasp.presets`.
        calc_defaults
            Default parameters for the recipe.
        calc_swaps
            Dictionary of custom kwargs for the Vasp calculator. Set a value to
            `None` to remove a pre-existing key entirely. For a list of available
            keys, refer to [quacc.calculators.vasp.vasp.Vasp][].
        report_mp_corrections
            Whether to report the Materials Project corrections in the results.
        additional_fields
            Additional fields to supply to the summarizer.
        copy_files
            Files to copy (and decompress) from source to the runtime directory.

        Returns
        -------
        None
        """
        self.input_atoms = input_atoms
        self.preset = preset
        self.calc_defaults = calc_defaults
        self.calc_swaps = calc_swaps
        self.report_mp_corrections = report_mp_corrections
        self.additional_fields = additional_fields
        self.copy_files = copy_files

        calc_flags = recursive_dict_merge(calc_defaults, calc_swaps)
        self.input_atoms.calc = Vasp(self.input_atoms, preset=preset, **calc_flags)

    def calculate(self) -> VaspSchema:
        """
        Base job function for VASP recipes.

        Returns
        -------
        VaspSchema
            Dictionary of results
        """
        final_atoms = run_calc(self.input_atoms, copy_files=self.copy_files)
        return vasp_summarize_run(
            final_atoms,
            report_mp_corrections=self.report_mp_corrections,
            additional_fields=self.additional_fields,
        )

    def optimize(
        self,
        opt_defaults: dict[str, Any] | None = None,
        opt_params: dict[str, Any] | None = None,
    ) -> VaspASESchema:
        """
        Base job function for VASP recipes with ASE optimizers.

        Parameters
        ----------
        opt_defaults
            Default arguments for the ASE optimizer.
        opt_params
            Dictionary of custom kwargs for [quacc.runners.ase.run_opt][]

        Returns
        -------
        VaspASESchema
            Dictionary of results
        """
        opt_flags = recursive_dict_merge(opt_defaults, opt_params)
        dyn = run_opt(self.input_atoms, copy_files=self.copy_files, **opt_flags)
        return summarize_vasp_opt_run(
            dyn,
            report_mp_corrections=self.report_mp_corrections,
            additional_fields=self.additional_fields,
        )
