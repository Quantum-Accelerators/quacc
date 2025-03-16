"""Schemas for storing ASE-based data."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from ase.io import read
from ase.vibrations.data import VibrationsData
from emmet.core.symmetry import PointGroupData
from pymatgen.io.ase import AseAtomsAdaptor

from quacc import QuaccDefault, __version__, get_settings
from quacc.atoms.core import get_final_atoms_from_dynamics
from quacc.schemas.atoms import atoms_to_metadata
from quacc.schemas.prep import prep_next_run
from quacc.schemas.thermo import ThermoSummarize
from quacc.utils.dicts import finalize_dict, recursive_dict_merge
from quacc.utils.files import get_uri

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Any, Literal

    from ase.atoms import Atoms
    from ase.io.trajectory import TrajectoryWriter
    from ase.md.md import MolecularDynamics
    from ase.optimize.optimize import Optimizer
    from ase.vibrations import Vibrations
    from maggma.core import Store

    from quacc.types import (
        DefaultSetting,
        DynSchema,
        OptSchema,
        RunSchema,
        VibSchema,
        VibThermoSchema,
    )


class Summarize:
    """
    Get tabulated results from an Atoms object and calculator and store them in a
    database-friendly format. This is meant to be compatible with all calculator types.
    """

    def __init__(
        self,
        directory: str | Path | None = None,
        move_magmoms: bool = False,
        additional_fields: dict[str, Any] | None = None,
    ) -> None:
        """
        Initialize the Summarize object.

        Parameters
        ----------
        directory
            Path to the directory where the calculation was run and results will be stored.
        move_magmoms
            Whether to move the final magmoms of the original Atoms object to the
            initial magmoms of the returned Atoms object, if relevant.
        additional_fields
            Additional fields to add to the task document.

        Returns
        -------
        None
        """
        self.directory = directory
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
        Get tabulated results from a standard ASE run.

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

        # Check and set up variables
        if not final_atoms.calc:
            msg = "ASE Atoms object has no attached calculator."
            raise ValueError(msg)
        if not final_atoms.calc.results:
            msg = "ASE Atoms object's calculator has no results."
            raise ValueError(msg)

        store = self._settings.STORE if store == QuaccDefault else store
        directory = self.directory or final_atoms.calc.directory

        # Generate input atoms metadata
        input_atoms_metadata = atoms_to_metadata(input_atoms) if input_atoms else {}

        # Generate the base of the task document
        inputs = {
            "parameters": final_atoms.calc.parameters,
            "nid": get_uri(directory).split(":")[0],
            "dir_name": directory,
            "input_atoms": input_atoms_metadata,
            "quacc_version": __version__,
        }
        results = {"results": final_atoms.calc.results}

        # Prepare atoms for the next run
        atoms_to_store = prep_next_run(final_atoms, move_magmoms=self.move_magmoms)

        # Generate final atoms metadata
        final_atoms_metadata = atoms_to_metadata(atoms_to_store) if final_atoms else {}

        # Create a dictionary of the inputs/outputs
        unsorted_task_doc = (
            final_atoms_metadata | inputs | results | self.additional_fields
        )

        return finalize_dict(
            unsorted_task_doc,
            directory=directory,
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
        Get tabulated results from an ASE optimization.

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

        # Check and set up variables
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
        directory = self.directory or final_atoms.calc.directory

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
        Get tabulated results from an ASE MD run.

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
        # Check and set up variables
        store = self._settings.STORE if store == QuaccDefault else store
        base_task_doc = self.opt(
            dyn, trajectory=trajectory, check_convergence=False, store=None
        )
        del base_task_doc["converged"]
        directory = self.directory or base_task_doc["dir_name"]

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
            directory=directory,
            gzip_file=self._settings.GZIP_FILES,
            store=store,
        )

    def neb(
        self,
        dyn: Optimizer,
        n_images: int,
        n_iter_return: int = -1,
        trajectory: TrajectoryWriter | list[Atoms] | None = None,
        store: Store | None | DefaultSetting = QuaccDefault,
    ) -> OptSchema:
        """
        Summarize the NEB run results and store them in a database-friendly format.

        Parameters
        ----------
        dyn
            ASE Optimizer object used for the NEB run.
        n_images
            Number of images in the NEB run.
        n_iter_return
            Number of iterations to return. If -1, all iterations are returned.
        trajectory
            Trajectory of the NEB run, either as a Trajectory object or a list of Atoms objects.
        store
            Maggma Store object to store the results in. Defaults to `QuaccSettings.STORE`.

        Returns
        -------
        OptSchema
            A dictionary containing the summarized NEB run results.
        """
        store = self._settings.STORE if store == QuaccDefault else store

        # Get trajectory
        if trajectory:
            atoms_trajectory = trajectory
        else:
            atoms_trajectory = read(dyn.trajectory.filename, index=":")  # type: ignore[union-attr]

        if n_iter_return == -1:
            atoms_trajectory = atoms_trajectory[-(n_images):]
        else:
            atoms_trajectory = _get_nth_iteration(
                atoms_trajectory,
                int(len(atoms_trajectory) / n_images),
                n_images,
                n_iter_return,
            )
        trajectory_results = [atoms.calc.results for atoms in atoms_trajectory]
        ts_index = (
            np.argmax(
                [
                    result["energy"]
                    for result in trajectory_results[-(n_images - 1) : -1]
                ]
            )
            + 1
        )
        ts_atoms = atoms_trajectory[ts_index]
        base_task_doc = atoms_to_metadata(atoms_trajectory[0])

        # Clean up the opt parameters
        parameters_opt = dyn.todict()
        parameters_opt.pop("logfile", None)
        parameters_opt.pop("restart", None)

        opt_fields = {
            "parameters_opt": parameters_opt,
            "trajectory": atoms_trajectory,
            "trajectory_results": trajectory_results,
            "ts_atoms": ts_atoms,
        }

        # Create a dictionary of the inputs/outputs
        unsorted_task_doc = base_task_doc | opt_fields | self.additional_fields

        return finalize_dict(
            unsorted_task_doc,
            directory=None,
            gzip_file=self._settings.GZIP_FILES,
            store=store,
        )


class VibSummarize:
    """
    Summarize an ASE Vibrations analysis.
    """

    def __init__(
        self,
        vib_object: Vibrations | VibrationsData,
        directory: str | Path | None = None,
        additional_fields: dict[str, Any] | None = None,
    ) -> None:
        """
        Initialize the Summarize object.

        Parameters
        ----------
        vib_object
            Instantiated ASE Vibrations object.
        directory
            Path to the directory where the results will be stored.
        additional_fields
            Additional fields to add to the task document.

        Returns
        -------
        None
        """
        self.vib_object = vib_object
        self.directory = directory
        self.additional_fields = additional_fields or {}
        self._settings = get_settings()

    def vib(
        self,
        is_molecule: bool = False,
        store: Store | None | DefaultSetting = QuaccDefault,
    ) -> VibSchema:
        """
        Get tabulated results from an ASE Vibrations object and store them in a database-
        friendly format.

        Parameters
        ----------
        is_molecule
            Whether the Atoms object is a molecule. If True, the vibrational modes are
            sorted by their absolute value and the 3N-5 or 3N-6 modes are taken. If False,
            all vibrational modes are taken.
        store
            Maggma Store object to store the results in. Defaults to `QuaccSettings.STORE`

        Returns
        -------
        VibSchema
            Dictionary representation of the task document
        """
        store = self._settings.STORE if store == QuaccDefault else store

        # Tabulate input parameters
        vib_freqs_raw = self.vib_object.get_frequencies().tolist()
        vib_energies_raw = self.vib_object.get_energies().tolist()
        if isinstance(self.vib_object, VibrationsData):
            atoms = self.vib_object._atoms
            directory = self.directory
            inputs = {"nid": get_uri(directory).split(":")[0], "dir_name": directory}
        else:
            atoms = self.vib_object.atoms
            directory = self.directory or atoms.calc.directory
            inputs = {
                "parameters": atoms.calc.parameters,
                "parameters_vib": {
                    "delta": self.vib_object.delta,
                    "direction": self.vib_object.direction,
                    "method": self.vib_object.method,
                    "ndof": self.vib_object.ndof,
                    "nfree": self.vib_object.nfree,
                },
                "nid": get_uri(directory).split(":")[0],
                "dir_name": directory,
                "quacc_version": __version__,
            }

        # Convert imaginary modes to negative values for DB storage
        for i, f in enumerate(vib_freqs_raw):
            if np.imag(f) > 0:
                vib_freqs_raw[i] = -np.abs(f)
                vib_energies_raw[i] = -np.abs(vib_energies_raw[i])
            else:
                vib_freqs_raw[i] = np.abs(f)
                vib_energies_raw[i] = np.abs(vib_energies_raw[i])

        # Get the true vibrational modes
        atoms_metadata = atoms_to_metadata(atoms)

        natoms = len(atoms)
        if natoms == 1:
            vib_freqs = []
            vib_energies = []
        elif is_molecule:
            is_linear = (
                PointGroupData()
                .from_molecule(AseAtomsAdaptor().get_molecule(atoms))
                .linear
                if atoms.pbc.any()
                else atoms_metadata["molecule_metadata"]["symmetry"]["linear"]
            )

            # Sort by absolute value
            vib_freqs_raw_sorted = vib_freqs_raw.copy()
            vib_energies_raw_sorted = vib_energies_raw.copy()
            vib_freqs_raw_sorted.sort(key=np.abs)
            vib_energies_raw_sorted.sort(key=np.abs)

            # Cut the 3N-5 or 3N-6 modes based on their absolute value
            n_modes = 3 * natoms - 5 if is_linear else 3 * natoms - 6
            vib_freqs = vib_freqs_raw_sorted[-n_modes:]
            vib_energies = vib_energies_raw_sorted[-n_modes:]
        else:
            vib_freqs = vib_freqs_raw
            vib_energies = vib_energies_raw

        imag_vib_freqs = [f for f in vib_freqs if f < 0]

        vib_results = {
            "results": {
                "imag_vib_freqs": imag_vib_freqs,
                "n_imag": len(imag_vib_freqs),
                "vib_energies": vib_energies,
                "vib_freqs": vib_freqs,
                "vib_energies_raw": vib_energies_raw,
                "vib_freqs_raw": vib_freqs_raw,
            }
        }
        unsorted_task_doc = (
            atoms_metadata | inputs | vib_results | self.additional_fields
        )

        return finalize_dict(
            unsorted_task_doc,
            directory=directory,
            gzip_file=self._settings.GZIP_FILES,
            store=store,
        )

    def vib_and_thermo(
        self,
        thermo_method: Literal["ideal_gas", "harmonic"],
        energy: float = 0.0,
        temperature: float = 298.15,
        pressure: float = 1.0,
        store: Store | None | DefaultSetting = QuaccDefault,
    ) -> VibThermoSchema:
        """
        Get tabulated results from an ASE Vibrations object and thermochemistry.

        Parameters
        ----------
        thermo_method
            Method to use for thermochemistry calculations. If None, no thermochemistry
            calculations are performed.
        energy
            Potential energy in eV used as the reference point for thermochemistry calculations.
        temperature
            Temperature in K for thermochemistry calculations.
        pressure
            Pressure in atm for thermochemistry calculations
        store
            Maggma Store object to store the results in. Defaults to `QuaccSettings.STORE`

        Returns
        -------
        VibThermoSchema
            Dictionary representation of the task document
        """
        store = self._settings.STORE if store == QuaccDefault else store

        atoms = (
            self.vib_object._atoms
            if isinstance(self.vib_object, VibrationsData)
            else self.vib_object.atoms
        )
        is_molecule = bool(thermo_method == "ideal_gas")

        # Generate vib data
        vib_schema = self.vib(is_molecule=is_molecule, store=None)
        directory = vib_schema["dir_name"]

        # Generate thermo data
        thermo_summary = ThermoSummarize(
            atoms,
            vib_schema["results"]["vib_freqs_raw"],
            energy=energy,
            directory=directory,
            additional_fields=self.additional_fields,
        )
        if thermo_method == "ideal_gas":
            thermo_schema = thermo_summary.ideal_gas(
                temperature=temperature, pressure=pressure, store=None
            )
        elif thermo_method == "harmonic":
            thermo_schema = thermo_summary.harmonic(
                temperature=temperature, pressure=pressure, store=None
            )
        else:
            raise ValueError(f"Unsupported thermo_method: {thermo_method}.")

        # Merge the vib and thermo data
        unsorted_task_doc = recursive_dict_merge(vib_schema, thermo_schema)

        return finalize_dict(
            unsorted_task_doc,
            directory=directory,
            gzip_file=self._settings.GZIP_FILES,
            store=store,
        )


def _get_nth_iteration(
    neb_trajectory: list[Atoms], n_iter: int, n_images: int, interval: int
) -> list[Atoms]:
    """
    Extract every nth iteration from the NEB trajectory.

    Parameters
    ----------
    neb_trajectory
        List of configurations (length: n_iter * n_images).
    n_iter
        Total number of iterations.
    n_images
        Number of images per iteration.
    interval
        Interval to get every nth iteration.

    Returns
    -------
    list[Atoms]
        List of configurations from every nth iteration.
    """
    result = []
    start_idx, end_idx = 0, 0
    for i in range(0, n_iter, interval):
        start_idx = i * n_images
        end_idx = start_idx + n_images

        result.extend(neb_trajectory[start_idx:end_idx])
    if end_idx < len(neb_trajectory) - 1:
        result.extend(neb_trajectory[-(n_images):])
    return result
