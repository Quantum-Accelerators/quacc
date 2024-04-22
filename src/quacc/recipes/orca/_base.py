"""Base jobs for ORCA."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.orca import ORCA, OrcaProfile, OrcaTemplate

from quacc import SETTINGS
from quacc.runners.ase import run_calc, run_opt
from quacc.schemas.cclib import cclib_summarize_run, summarize_cclib_opt_run
from quacc.utils.dicts import recursive_dict_merge
from quacc.utils.lists import merge_list_params

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.cclib import cclibASEOptSchema, cclibSchema
    from quacc.utils.files import Filenames, SourceDirectory

_LABEL = OrcaTemplate()._label  # skipcq: PYL-W0212
LOG_FILE = f"{_LABEL}.out"
GEOM_FILE = f"{_LABEL}.xyz"


class RunAndSummarize:
    """Handler for running and summarizing ORCA calculations."""

    def __init__(
        self,
        atoms: Atoms,
        charge: int,
        spin_multiplicity: int,
        default_inputs: list[str] | None = None,
        default_blocks: list[str] | None = None,
        input_swaps: list[str] | None = None,
        block_swaps: list[str] | None = None,
        additional_fields: dict[str, Any] | None = None,
        copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
        **calc_kwargs,
    ) -> None:
        """
        Initialize the class.

        Parameters
        ----------
        atoms
            Atoms object
        charge
            Charge of the system.
        spin_multiplicity
            Multiplicity of the system.
        default_inputs
            Default input parameters.
        default_blocks
            Default block input parameters.
        input_swaps
            List of orcasimpleinput swaps for the calculator. To remove entries
            from the defaults, put a `#` in front of the name.
        block_swaps
            List of orcablock swaps for the calculator. To remove entries
            from the defaults, put a `#` in front of the name.
        additional_fields
            Any additional fields to supply to the summarizer.
        copy_files
            Files to copy (and decompress) from source to the runtime directory.
        **calc_kwargs
            Any other keyword arguments to pass to the `ORCA` calculator.
        """
        self.atoms = atoms
        self.charge = charge
        self.spin_multiplicity = spin_multiplicity
        self.default_inputs = default_inputs
        self.default_blocks = default_blocks
        self.input_swaps = input_swaps
        self.block_swaps = block_swaps
        self.additional_fields = additional_fields
        self.copy_files = copy_files
        self.calc_kwargs = calc_kwargs
        self.atoms.calc = self._prep_calculator()

    def calculate(self) -> cclibSchema:
        """
        Base job function for ORCA recipes.

        Returns
        -------
        cclibSchema
            Dictionary of results
        """
        atoms = run_calc(self.atoms, geom_file=GEOM_FILE, copy_files=self.copy_files)
        return cclib_summarize_run(
            atoms, LOG_FILE, additional_fields=self.additional_fields
        )

    def optimize(
        self,
        opt_defaults: dict[str, Any] | None = None,
        opt_params: dict[str, Any] | None = None,
    ) -> cclibASEOptSchema:
        """
        Base job function for ORCA recipes with ASE optimizer.

        Parameters
        ----------
        opt_defaults
            Default arguments for the ASE optimizer.
        opt_params
            Dictionary of custom kwargs for [quacc.runners.ase.run_opt][]

        Returns
        -------
        cclibASEOptSchema
            Dictionary of results
        """
        opt_flags = recursive_dict_merge(opt_defaults, opt_params)
        dyn = run_opt(self.atoms, copy_files=self.copy_files, **opt_flags)
        return summarize_cclib_opt_run(
            dyn, LOG_FILE, additional_fields=self.additional_fields
        )

    def _prep_calculator(self) -> ORCA:
        """
        Prepare the ORCA calculator.

        Returns
        -------
        ORCA
            The ORCA calculator
        """
        inputs = merge_list_params(self.default_inputs, self.input_swaps)
        blocks = merge_list_params(self.default_blocks, self.block_swaps)
        if "xyzfile" not in inputs:
            inputs.append("xyzfile")
        orcasimpleinput = " ".join(inputs)
        orcablocks = "\n".join(blocks)

        return ORCA(
            profile=OrcaProfile(SETTINGS.ORCA_CMD),
            charge=self.charge,
            mult=self.spin_multiplicity,
            orcasimpleinput=orcasimpleinput,
            orcablocks=orcablocks,
            **self.calc_kwargs,
        )
