"""Schemas for storing ASE-based data."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.io import read

from quacc import QuaccDefault, __version__, get_settings
from quacc.atoms.core import get_final_atoms_from_dynamics
from quacc.schemas.atoms import atoms_to_metadata
from quacc.schemas.prep import prep_next_run
from quacc.utils.dicts import finalize_dict
from quacc.utils.files import get_uri

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms
    from ase.md.md import MolecularDynamics
    from ase.optimize.optimize import Optimizer
    from maggma.core import Store

    from quacc.types import DefaultSetting, DynSchema, OptSchema, RunSchema


class Summarize:
    def __init__(
        self,
        charge_and_multiplicity: tuple[int, int] | None = None,
        move_magmoms: bool = False,
        additional_fields: dict[str, Any] | None = None,
    ) -> None:
        """
        Get tabulated results from an Atoms object and calculator and store them in a
        database-friendly format. This is meant to be compatible with all calculator types.

        Parameters
        ----------
        charge_and_multiplicity
            Charge and spin multiplicity of the Atoms object, only used for Molecule
            metadata.
        move_magmoms
            Whether to move the final magmoms of the original Atoms object to the
            initial magmoms of the returned Atoms object.
        additional_fields
            Additional fields to add to the task document.
        """
        self.charge_and_multiplicity = charge_and_multiplicity
        self.move_magmoms = move_magmoms
        self.additional_fields = additional_fields or {}
        self._settings = get_settings()

    def run(
        self,
        final_atoms: Atoms,
        input_atoms: Atoms,
        store: Store | None | DefaultSetting = QuaccDefault,
    ) -> RunSchema:
        """
        Get tabulated results from an Atoms object and calculator and store them in a
        database-friendly format. This is meant to be compatible with all calculator types.

        Parameters
        ----------
        final_atoms
            ASE Atoms following a calculation. A calculator must be attached.
        input_atoms
            Input ASE Atoms object to store.
        store
            Maggma Store object to store the results in. Defaults to `QuaccSettings.STORE`

        Returns
        -------
        RunSchema
            Dictionary representation of the task document
        """
        if not final_atoms.calc:
            msg = "ASE Atoms object has no attached calculator."
            raise ValueError(msg)
        if not final_atoms.calc.results:
            msg = "ASE Atoms object's calculator has no results."
            raise ValueError(msg)
        store = self._settings.STORE if store == QuaccDefault else store

        directory = final_atoms.calc.directory
        uri = get_uri(directory)

        if input_atoms:
            input_atoms_metadata = atoms_to_metadata(
                input_atoms,
                charge_and_multiplicity=self.charge_and_multiplicity,
                store_pmg=False,
            )
        else:
            input_atoms_metadata = {}

        inputs = {
            "parameters": final_atoms.calc.parameters,
            "nid": uri.split(":")[0],
            "dir_name": directory,
            "input_atoms": input_atoms_metadata,
            "quacc_version": __version__,
        }
        results = {"results": final_atoms.calc.results}

        atoms_to_store = prep_next_run(final_atoms, move_magmoms=self.move_magmoms)

        if final_atoms:
            final_atoms_metadata = atoms_to_metadata(
                atoms_to_store, charge_and_multiplicity=self.charge_and_multiplicity
            )
        else:
            final_atoms_metadata = {}

        unsorted_task_doc = (
            final_atoms_metadata | inputs | results | self.additional_fields
        )

        return finalize_dict(
            unsorted_task_doc,
            directory,
            gzip_file=self._settings.GZIP_FILES,
            store=store,
        )

    def opt(
        self,
        dyn: Optimizer,
        trajectory: list[Atoms] | None = None,
        check_convergence: bool | DefaultSetting = QuaccDefault,
        store: Store | None | DefaultSetting = QuaccDefault,
    ) -> OptSchema:
        """
        Get tabulated results from an ASE Atoms trajectory and store them in a database-
        friendly format. This is meant to be compatible with all calculator types.

        Parameters
        ----------
        dyn
            ASE Optimizer object.
        trajectory
            ASE Trajectory object or list[Atoms] from reading a trajectory file. If
            None, the trajectory must be found in `dyn.trajectory.filename`.
        check_convergence
            Whether to check the convergence of the calculation. Defaults to True in
            settings.
        store
            Maggma Store object to store the results in. Defaults to `QuaccSettings.STORE`

        Returns
        -------
        OptSchema
            Dictionary representation of the task document
        """
        check_convergence = (
            self._settings.CHECK_CONVERGENCE
            if check_convergence == QuaccDefault
            else check_convergence
        )
        store = self._settings.STORE if store == QuaccDefault else store

        # Get trajectory
        if trajectory:
            atoms_trajectory = trajectory
        else:
            atoms_trajectory = read(dyn.trajectory.filename, index=":")  # type: ignore[union-attr]

        trajectory_results = [atoms.calc.results for atoms in atoms_trajectory]

        initial_atoms = atoms_trajectory[0]
        final_atoms = get_final_atoms_from_dynamics(dyn)
        directory = final_atoms.calc.directory

        # Check convergence
        is_converged = dyn.converged()
        if check_convergence and not is_converged:
            msg = f"Optimization did not converge. Refer to {directory}"
            raise RuntimeError(msg)

        # Base task doc
        base_task_doc = self.run(final_atoms, initial_atoms, store=None)

        # Clean up the opt parameters
        parameters_opt = dyn.todict()
        parameters_opt.pop("logfile", None)
        parameters_opt.pop("restart", None)

        opt_fields = {
            "parameters_opt": parameters_opt,
            "converged": is_converged,
            "trajectory": atoms_trajectory,
            "trajectory_results": trajectory_results,
        }

        # Create a dictionary of the inputs/outputs
        unsorted_task_doc = base_task_doc | opt_fields | self.additional_fields

        return finalize_dict(
            unsorted_task_doc,
            directory,
            gzip_file=self._settings.GZIP_FILES,
            store=store,
        )

    def md(
        self,
        dyn: MolecularDynamics,
        trajectory: list[Atoms] | None = None,
        store: Store | None | DefaultSetting = QuaccDefault,
    ) -> DynSchema:
        """
        Get tabulated results from an ASE Atoms trajectory and store them in a database-
        friendly format. This is meant to be compatible with all calculator types.

        Parameters
        ----------
        dyn
            ASE MolecularDynamics object.
        trajectory
            ASE Trajectory object or list[Atoms] from reading a trajectory file. If
            None, the trajectory must be found in `dyn.trajectory.filename`.
        store
            Maggma Store object to store the results in. Defaults to `QuaccSettings.STORE`

        Returns
        -------
        DynSchema
            Dictionary representation of the task document
        """
        store = self._settings.STORE if store == QuaccDefault else store
        base_task_doc = self.opt(
            dyn, trajectory=trajectory, check_convergence=False, store=None
        )
        del base_task_doc["converged"]

        # Clean up the opt parameters
        parameters_md = base_task_doc.pop("parameters_opt")
        parameters_md.pop("logfile", None)

        trajectory_log = []
        for t, atoms in enumerate(base_task_doc["trajectory"]):
            trajectory_log.append(
                {
                    "kinetic_energy": atoms.get_kinetic_energy(),
                    "temperature": atoms.get_temperature(),
                    "time": t * parameters_md["timestep"],
                }
            )

        md_fields = {"parameters_md": parameters_md, "trajectory_log": trajectory_log}

        # Create a dictionary of the inputs/outputs
        unsorted_task_doc = base_task_doc | md_fields | self.additional_fields

        return finalize_dict(
            unsorted_task_doc,
            base_task_doc["dir_name"],
            gzip_file=self._settings.GZIP_FILES,
            store=store,
        )
