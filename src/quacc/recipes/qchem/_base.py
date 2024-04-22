"""Base jobs for Q-Chem."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc.calculators.qchem import QChem
from quacc.runners.ase import run_calc, run_opt
from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import OptSchema, RunSchema
    from quacc.utils.files import Filenames, SourceDirectory


class RunAndSummarize:
    """Run and summarize a Q-Chem calculation."""

    def __init__(
        self,
        input_atoms: Atoms,
        charge: int = 0,
        spin_multiplicity: int = 1,
        calc_defaults: dict[str, Any] | None = None,
        calc_swaps: dict[str, Any] | None = None,
        additional_fields: dict[str, Any] | None = None,
        copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
        calc_defaults
            The default parameters for the recipe.
        calc_swaps
            Dictionary of custom kwargs for the Q-Chem calculator. Set a value to `quacc.Remove` to
            remove a pre-existing key entirely. For a list of available keys, refer to the
            `quacc.calculators.qchem.qchem.QChem` calculator.
        additional_fields
            Any additional fields to set in the summary.
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

        calc_flags = recursive_dict_merge(self.calc_defaults, self.calc_swaps)
        self.input_atoms.calc = QChem(
            self.input_atoms,
            charge=self.charge,
            spin_multiplicity=self.spin_multiplicity,
            **calc_flags,
        )

    def calculate(self) -> RunSchema:
        """
        Base job function used for Q-Chem recipes that don't rely on ASE optimizers or other
        ASE dynamics classes.

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

    def optimize(
        self,
        opt_defaults: dict[str, Any] | None = None,
        opt_params: dict[str, Any] | None = None,
    ) -> OptSchema:
        """
        Base function for Q-Chem recipes that involve ASE optimizers.

        Parameters
        ----------
        opt_defaults
            Default arguments for the ASE optimizer.
        opt_params
            Dictionary of custom kwargs for [quacc.runners.ase.run_opt][]

        Returns
        -------
        OptSchema
            Dictionary of results from [quacc.schemas.ase.summarize_opt_run][]
        """
        # TODO:
        #   - passing initial Hessian?

        opt_flags = recursive_dict_merge(opt_defaults, opt_params)
        dyn = run_opt(self.input_atoms, copy_files=self.copy_files, **opt_flags)

        return summarize_opt_run(
            dyn,
            charge_and_multiplicity=(self.charge, self.spin_multiplicity),
            additional_fields=self.additional_fields,
        )
