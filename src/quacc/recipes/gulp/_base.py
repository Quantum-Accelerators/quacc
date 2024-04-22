"""Base jobs for GULP."""

from __future__ import annotations

import logging
import os
from typing import TYPE_CHECKING

from ase.calculators.gulp import GULP

from quacc import SETTINGS
from quacc.runners.ase import run_calc
from quacc.schemas.ase import summarize_run
from quacc.utils.lists import merge_list_params

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import RunSchema
    from quacc.utils.files import Filenames, SourceDirectory
logger = logging.getLogger(__name__)

GEOM_FILE_PBC = "gulp.cif"
GEOM_FILE_NOPBC = "gulp.xyz"
GULP_CMD = f"{SETTINGS.GULP_CMD} < gulp.gin > gulp.got"


class RunAndSummarize:
    """Run and summarize a GULP run."""

    def __init__(
        self,
        input_atoms: Atoms,
        library: str | None = None,
        keyword_defaults: list[str] | None = None,
        option_defaults: list[str] | None = None,
        keyword_swaps: list[str] | None = None,
        option_swaps: list[str] | None = None,
        additional_fields: dict[str, Any] | None = None,
        copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
        **calc_kwargs,
    ) -> None:
        """
        Initialize the class.

        Parameters
        ----------
        input_atoms
            Atoms object
        library
            Filename of the potential library file, if required.
        keyword_defaults
            Default `keywords` for calculator.
        option_defaults
            Default `options` for calculator.
        keyword_swaps
            List of custom `keyword` kwargs for the GULP calculator. To remove entries
            from the defaults, put a `#` in front of the name. For a list of
            available keys, refer to the `ase.calculators.gulp.GULP` calculator.
        option_swaps
            Dictionary of custom `options` kwargs for the GULP calculator. To remove entries
            from the defaults, put a `#` in front of the name. For a list of
            available keys, refer to the `ase.calculators.gulp.GULP` calculator.
        additional_fields
            Additional field to supply to the summarizer.
        copy_files
            Files to copy (and decompress) from source to the runtime directory.
        **calc_kwargs
            Any other keyword arguments to pass to the `GULP` calculator.

        Returns
        -------
        None
        """
        self.input_atoms = input_atoms
        self.library = library
        self.keyword_defaults = keyword_defaults
        self.option_defaults = option_defaults
        self.keyword_swaps = keyword_swaps
        self.option_swaps = option_swaps
        self.additional_fields = additional_fields
        self.copy_files = copy_files
        self.calc_kwargs = calc_kwargs
        self._prep_calculator()

    def calculate(self) -> RunSchema:
        """
        Base job function for GULP recipes.

        Returns
        -------
        RunSchema
            Dictionary of results from [quacc.schemas.ase.summarize_run][]
        """
        final_atoms = run_calc(
            self.input_atoms,
            geom_file=GEOM_FILE_PBC if self.input_atoms.pbc.any() else GEOM_FILE_NOPBC,
            copy_files=self.copy_files,
        )

        if (
            SETTINGS.CHECK_CONVERGENCE
            and "opti" in self.gulp_keywords
            and not final_atoms.calc.get_opt_state()
        ):
            msg = "Optimization did not converge."
            raise RuntimeError(msg)

        return summarize_run(
            final_atoms, self.input_atoms, additional_fields=self.additional_fields
        )

    def _prep_calculator(self) -> None:
        """
        Prepare the calculator for the GULP calculation.

        Returns
        -------
        None
        """
        self.keyword_defaults = self.keyword_defaults or []

        if not self.input_atoms.pbc.any():
            if "opti" in self.keyword_defaults and "conv" not in self.keyword_defaults:
                self.keyword_defaults += ["conv"]
            self.keyword_defaults = [
                k for k in self.keyword_defaults if k not in ["gwolf", "conp"]
            ]

        self.option_defaults += [
            (
                f"output cif {GEOM_FILE_PBC}"
                if self.input_atoms.pbc.any()
                else f"output xyz {GEOM_FILE_NOPBC}"
            )
        ]

        keywords = merge_list_params(self.keyword_defaults, self.keyword_swaps)
        options = merge_list_params(self.option_defaults, self.option_swaps)

        self.gulp_keywords = " ".join(keywords)
        self.gulp_options = list(options)

        if SETTINGS.GULP_LIB:
            os.environ["GULP_LIB"] = str(SETTINGS.GULP_LIB)
        self.input_atoms.calc = GULP(
            command=GULP_CMD,
            keywords=self.gulp_keywords,
            options=self.gulp_options,
            library=self.library,
            **self.calc_kwargs,
        )
