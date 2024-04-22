"""Base jobs for Onetep."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.onetep import Onetep, OnetepProfile

from quacc import SETTINGS
from quacc.runners.ase import run_calc, run_opt
from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import RunSchema
    from quacc.utils.files import Filenames, SourceDirectory


class RunAndSummarize:
    """Run and summarize a ONETEP calculation."""

    def __init__(
        self,
        input_atoms: Atoms,
        calc_defaults: dict[str, Any] | None = None,
        calc_swaps: dict[str, Any] | None = None,
        additional_fields: dict[str, Any] | None = None,
        copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    ) -> None:
        """
        Initialize the class.

        Parameters
        ----------
        input_atoms
            Atoms object
        calc_defaults
            The default calculator parameters.
        calc_swaps
            Custom kwargs for the ONETEP calculator. Set a value to
            `quacc.Remove` to remove a pre-existing key entirely. For a list of available
            keys, refer to the [ase.calculators.onetep.Onetep][] calculator.
        additional_fields
            Any additional fields to supply to the summarizer.
        copy_files
            Files to copy (and decompress) from source to the runtime directory.

        Returns
        -------
        None
        """
        self.input_atoms = input_atoms
        self.calc_defaults = calc_defaults
        self.calc_swaps = calc_swaps
        self.additional_fields = additional_fields
        self.copy_files = copy_files

        calc_flags = recursive_dict_merge(calc_defaults, calc_swaps)
        self.input_atoms.calc = Onetep(
            pseudo_path=str(SETTINGS.ONETEP_PP_PATH)
            if SETTINGS.ONETEP_PP_PATH
            else ".",
            parallel_info=SETTINGS.ONETEP_PARALLEL_CMD,
            profile=OnetepProfile(
                SETTINGS.ONETEP_CMD
            ),  # TODO: If the ASE merge is successful, we need to change ONETEP_PARALLEL_CMD to a list[str] and remove parallel info.
            # If we also have access to post_args we can point not to the binary but to the launcher which takes -t nthreads as a post_args
            **calc_flags,
        )

    def calculate(self) -> RunSchema:
        """
        Base function to carry out Onetep recipes.

        Returns
        -------
        RunSchema
            Dictionary of results from [quacc.schemas.ase.summarize_run][]
        """
        final_atoms = run_calc(self.input_atoms, copy_files=self.copy_files)
        return summarize_run(
            final_atoms, self.input_atoms, additional_fields=self.additional_fields
        )

    def run_and_summarize_opt(
        self,
        opt_defaults: dict[str, Any] | None = None,
        opt_params: dict[str, Any] | None = None,
    ) -> RunSchema:
        """
        Base function to carry out Onetep recipes with ASE optimizers.

        Parameters
        ----------
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
        dyn = run_opt(self.input_atoms, copy_files=self.copy_files, **opt_flags)
        return summarize_opt_run(dyn, additional_fields=self.additional_fields)
