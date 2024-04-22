"""Base jobs for espresso."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from ase.atoms import Atoms
from ase.io.espresso import Namelist
from ase.io.espresso_namelist.keys import ALL_KEYS

from quacc.calculators.espresso.espresso import (
    Espresso,
    EspressoProfile,
    EspressoTemplate,
)
from quacc.calculators.espresso.utils import (
    prepare_copy_files,
    remove_conflicting_kpts_kspacing,
)
from quacc.runners.ase import run_calc, run_opt
from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from quacc.schemas._aliases.ase import RunSchema
    from quacc.utils.files import Filenames, SourceDirectory


class RunAndSummarize:
    """Run and summarize a Quantum Espresso calculation."""

    def __init__(
        self,
        input_atoms: Atoms | None = None,
        preset: str | None = None,
        template: EspressoTemplate | None = None,
        profile: EspressoProfile | None = None,
        calc_defaults: dict[str, Any] | None = None,
        calc_swaps: dict[str, Any] | None = None,
        parallel_info: dict[str, Any] | None = None,
        additional_fields: dict[str, Any] | None = None,
        copy_files: (
            SourceDirectory
            | list[SourceDirectory]
            | dict[SourceDirectory, Filenames]
            | None
        ) = None,
    ) -> None:
        """
        Initialize the class.

        Parameters
        ----------
        atoms
            Atoms object
        preset
            Name of the preset to use
        template
            EspressoTemplate to use
        profile
            EspressoProfile to use
        calc_defaults
            The default calculator parameters.
        calc_swaps
            Custom kwargs for the espresso calculator. Set a value to
            `quacc.Remove` to remove a pre-existing key entirely. For a list of available
            keys, refer to the [ase.calculators.espresso.Espresso][] calculator.
        parallel_info
            Dictionary of parallelization information.
        additional_fields
            Any additional fields to supply to the summarizer.
        copy_files
            Files to copy (and decompress) from source to the runtime directory.

        Returns
        -------
        None
        """
        self.input_atoms = input_atoms
        self.preset = preset
        self.template = template
        self.profile = profile
        self.calc_defaults = calc_defaults
        self.calc_swaps = calc_swaps
        self.parallel_info = parallel_info
        self.additional_fields = additional_fields
        self.copy_files = copy_files

        self._prepare_calc()
        self._prepare_copy()

    def calculate(self) -> RunSchema:
        """
        Base function to carry out espresso recipes.

        Parameters
        ----------

        Returns
        -------
        RunSchema
            Dictionary of results from [quacc.schemas.ase.summarize_run][]
        """
        geom_file = self.template.outputname if self.template.binary == "pw" else None
        final_atoms = run_calc(
            self.input_atoms, geom_file=geom_file, copy_files=self.copy_files
        )
        return summarize_run(
            final_atoms,
            self.input_atoms,
            move_magmoms=True,
            additional_fields=self.additional_fields,
        )

    def optimize(
        self,
        relax_cell: bool = False,
        opt_defaults: dict[str, Any] | None = None,
        opt_params: dict[str, Any] | None = None,
    ) -> RunSchema:
        """
        Base function to carry out espresso recipes with ASE optimizers.

        Parameters
        ----------
        relax_cell
            Whether to relax the cell or not.
        opt_defaults
            The default optimization parameters.
        opt_params
            Dictionary of parameters to pass to the optimizer. pass "optimizer"
            to change the optimizer being used. "fmax" and "max_steps" are commonly
            used keywords. See the ASE documentation for more information.

        Returns
        -------
        RunSchema
            Dictionary of results from [quacc.schemas.ase.summarize_run][]
        """
        opt_flags = recursive_dict_merge(opt_defaults, opt_params)
        dyn = run_opt(
            self.input_atoms,
            relax_cell=relax_cell,
            copy_files=self.updated_copy_files,
            **opt_flags,
        )
        return summarize_opt_run(
            dyn, move_magmoms=True, additional_fields=self.additional_fields
        )

    def _prepare_calc(self) -> None:
        """
        Commonly used preparation function to merge parameters
        and attach an Espresso calculator accordingly.

        -------
        None
        """
        self.input_atoms = Atoms() if self.input_atoms is None else self.input_atoms
        self.calc_defaults = self.calc_defaults or {}
        self.calc_swaps = self.calc_swaps or {}

        self.calc_defaults["input_data"] = Namelist(
            self.calc_defaults.get("input_data")
        )
        self.calc_swaps["input_data"] = Namelist(self.calc_swaps.get("input_data"))

        binary = self.template.binary if self.template else "pw"

        if binary in ALL_KEYS:
            self.calc_defaults["input_data"].to_nested(
                binary=binary, **self.calc_defaults
            )
            self.calc_swaps["input_data"].to_nested(binary=binary, **self.calc_swaps)

        self.calc_defaults = remove_conflicting_kpts_kspacing(
            self.calc_defaults, self.calc_swaps
        )

        self.calc_flags = recursive_dict_merge(self.calc_defaults, self.calc_swaps)

        self.input_atoms.calc = Espresso(
            input_atoms=self.input_atoms,
            preset=self.preset,
            parallel_info=self.parallel_info,
            template=self.template,
            profile=self.profile,
            **self.calc_flags,
        )

    def _prepare_copy(self) -> None:
        """
        Function that will prepare the files to copy.

        Returns
        -------
        None
        """
        if isinstance(self.copy_files, (str, Path)):
            self.copy_files = [self.copy_files]

        if isinstance(self.copy_files, list):
            exact_files_to_copy = prepare_copy_files(
                self.calc_flags, binary=self.input_atoms.calc.template.binary
            )
            self.copy_files = {
                source: exact_files_to_copy for source in self.copy_files
            }
