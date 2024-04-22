"""Base jobs for Psi4."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.psi4 import Psi4

from quacc.runners.ase import run_calc
from quacc.schemas.ase import summarize_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import RunSchema
    from quacc.utils.files import Filenames, SourceDirectory


class RunAndSummarize:
    """Run and summarize a PSI4 calculation."""

    def __init__(
        self,
        input_atoms: Atoms,
        charge: int,
        spin_multiplicity: int,
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
        charge
            Charge of the system.
        spin_multiplicity
            Multiplicity of the system.
        calc_defaults
            The default calculator parameters.
        calc_swaps
            Custom kwargs for the Psi4 calculator. Set a value to
            `quacc.Remove` to remove a pre-existing key entirely. For a list of available
            keys, refer to the [ase.calculators.psi4.Psi4][] calculator.
        additional_fields
            Any additional fields to supply to the summarizer.
        copy_files
            Files to copy (and decompress) from source to the runtime directory.

        Returns
        -------
        None
        """
        self.input_atoms = input_atoms
        self.charge = charge
        self.spin_multiplicity = spin_multiplicity
        self.calc_defaults = calc_defaults
        self.calc_swaps = calc_swaps
        self.additional_fields = additional_fields
        self.copy_files = copy_files
        self._prepare_calc()

    def calculate(self) -> RunSchema:
        """
        Base function to carry out Psi4 recipes.

        Returns
        -------
        RunSchema
            Dictionary of results from [quacc.schemas.ase.summarize_run][]
        """
        final_atoms = run_calc(self.input_atoms, copy_files=self.copy_files)
        return summarize_run(
            final_atoms,
            self.input_atoms,
            charge_and_multiplicity=(self.charge, self.spin_multiplicity),
            additional_fields=self.additional_fields,
        )

    def _prepare_calc(self) -> None:
        """
        Prepare the PSI4 calculator.

        Returns
        -------
        None
        """

        calc_flags = recursive_dict_merge(self.calc_defaults, self.calc_swaps)
        self.input_atoms.calc = Psi4(**calc_flags)
