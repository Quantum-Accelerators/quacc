"""
Utility functions for running ASE calculators
"""
from __future__ import annotations

import os
import warnings
from tempfile import mkdtemp

import numpy as np
from ase import Atoms
from ase.io import read, trajectory
from ase.optimize import (
    BFGS,
    FIRE,
    LBFGS,
    BFGSLineSearch,
    GPMin,
    LBFGSLineSearch,
    MDMin,
)
from ase.thermochemistry import IdealGasThermo
from ase.units import invcm
from ase.vibrations import Vibrations
from monty.io import zopen
from monty.os.path import zpath
from monty.shutil import copy_r, gzip_dir
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath

from quacc import SETTINGS
from quacc.schemas.atoms import atoms_to_metadata
from quacc.util.atoms import copy_atoms
from quacc.util.files import copy_decompress


def run_calc(
    atoms: Atoms,
    geom_file: str | None = None,
    scratch_dir: str = SETTINGS.SCRATCH_DIR,
    gzip: bool = SETTINGS.GZIP_FILES,
    copy_files: list[str] = None,
) -> Atoms:
    """
    Run a calculation in a scratch directory and copy the results back to the
    original directory. This can be useful if file I/O is slow in the working
    directory, so long as file transfer speeds are reasonable.

    This is a wrapper around atoms.get_potential_energy(). Note: This
    function does not modify the atoms object in-place.

    Parameters
    ----------
    atoms
        The Atoms object to run the calculation on.
    geom_file
        The filename of the log file that contains the output geometry, used
        to update the atoms object's positions and cell after a job. It is better
        to specify this rather than relying on ASE's atoms.get_potential_energy()
        function to update the positions, as this varies between codes.
    scratch_dir
        Path where a tmpdir should be made for running the calculation. If None,
        the current working directory will be used.
    gzip
        Whether to gzip the output files.
    copy_files
        Filenames to copy from source to scratch directory.

    Returns
    -------
    .Atoms
        The updated .Atoms object,
    """

    if atoms.calc is None:
        raise ValueError("Atoms object must have attached calculator.")
    atoms = copy_atoms(atoms)
    cwd = os.getcwd()
    scratch_dir = scratch_dir or cwd
    symlink = os.path.join(cwd, "tmp_dir")

    if not os.path.exists(scratch_dir):
        os.makedirs(scratch_dir)

    tmpdir = mkdtemp(prefix="quacc-tmp", dir=scratch_dir)

    if os.name != "nt":
        if os.path.islink(symlink):
            os.unlink(symlink)
        os.symlink(tmpdir, symlink)

    # Copy files to scratch and decompress them if needed
    if copy_files:
        copy_decompress(copy_files, tmpdir)

    # Run calculation via get_potential_energy()
    os.chdir(tmpdir)
    atoms.get_potential_energy()
    os.chdir(cwd)

    # Gzip files in tmpdir
    if gzip:
        gzip_dir(tmpdir)

    # Copy files back to run_dir
    copy_r(tmpdir, cwd)

    # Remove symlink
    if os.path.islink(symlink):
        os.remove(symlink)

    # Some ASE calculators do not update the atoms object in-place with
    # a call to .get_potential_energy(). This is a workaround to ensure
    # that the atoms object is updated with the correct positions, cell,
    # and magmoms.
    if geom_file and os.path.exists(zpath(geom_file)):
        # Note: We have to be careful to make sure we don't lose the
        # converged magnetic moments, if present. That's why we simply
        # update the positions and cell in-place.
        atoms_new = read(zpath(geom_file))
        if isinstance(atoms_new, list):
            atoms_new = atoms_new[-1]

        # Make sure the atom indices didn't get updated somehow (sanity check)
        if (
            np.array_equal(atoms_new.get_atomic_numbers(), atoms.get_atomic_numbers())
            is False
        ):
            raise ValueError("Atomic numbers do not match between atoms and geom_file.")

        atoms.positions = atoms_new.positions
        atoms.cell = atoms_new.cell

    return atoms


def run_ase_opt(
    atoms: Atoms,
    fmax: float = 0.01,
    max_steps: int = 500,
    optimizer: str = "FIRE",
    opt_kwargs: dict = None,
    scratch_dir: str = SETTINGS.SCRATCH_DIR,
    gzip: bool = SETTINGS.GZIP_FILES,
    copy_files: list[str] = None,
) -> trajectory:
    """
    Run an ASE-based optimization in a scratch directory and copy the results
    back to the original directory. This can be useful if file I/O is slow in
    the working directory, so long as file transfer speeds are reasonable.

    This is a wrapper around the optimizers in ASE. Note: This function does
    not modify the atoms object in-place.

    Parameters
    ----------
    atoms
        The Atoms object to run the calculation on.
    fmax
        Tolerance for the force convergence (in eV/A).
    max_steps
        Maximum number of steps to take.
    optimizer
        Name of optimizer class to use.
    opt_kwargs
        Dictionary of kwargs for the optimizer.
    scratch_dir
        Path where a tmpdir should be made for running the calculation. If None,
        the current working directory will be used.
    gzip
        Whether to gzip the output files.
    copy_files
        Filenames to copy from source to scratch directory.

    Returns
    -------
    traj
        The ASE trajectory object.
    """

    if atoms.calc is None:
        raise ValueError("Atoms object must have attached calculator.")

    atoms = copy_atoms(atoms)
    cwd = os.getcwd()
    scratch_dir = scratch_dir or cwd
    symlink = os.path.join(cwd, "tmp_dir")
    opt_kwargs = opt_kwargs or {}

    if not os.path.exists(scratch_dir):
        os.makedirs(scratch_dir)

    opt_kwargs["trajectory"] = "opt.traj"
    opt_kwargs["restart"] = "opt.pckl"

    # Get optimizer
    if optimizer.lower() == "bfgs":
        opt_class = BFGS
    elif optimizer.lower() == "bfgslinesearch":
        opt_class = BFGSLineSearch
    elif optimizer.lower() == "lbfgs":
        opt_class = LBFGS
    elif optimizer.lower() == "lbfgslinesearch":
        opt_class = LBFGSLineSearch
    elif optimizer.lower() == "gpmin":
        opt_class = GPMin
    elif optimizer.lower() == "mdmin":
        opt_class = MDMin
    elif optimizer.lower() == "fire":
        opt_class = FIRE
    else:
        raise ValueError(f"Unknown optimizer: {optimizer}")

    tmpdir = mkdtemp(prefix="quacc-tmp", dir=scratch_dir)

    if os.name != "nt":
        if os.path.islink(symlink):
            os.unlink(symlink)
        os.symlink(tmpdir, symlink)

    # Copy files to scratch and decompress them if needed
    if copy_files:
        copy_decompress(copy_files, tmpdir)

    # Run calculation
    os.chdir(tmpdir)
    dyn = opt_class(atoms, **opt_kwargs)
    dyn.run(fmax=fmax, steps=max_steps)
    os.chdir(cwd)

    # Check convergence
    if not dyn.converged:
        raise ValueError("Optimization did not converge.")

    # Read trajectory
    traj = read(os.path.join(tmpdir, "opt.traj"), index=":")

    # Gzip files in tmpdir
    if gzip:
        gzip_dir(tmpdir)

    # Copy files back to run_dir
    copy_r(tmpdir, cwd)

    # Remove symlink
    if os.path.islink(symlink):
        os.remove(symlink)

    os.chdir(cwd)

    return traj


def run_ase_vib(
    atoms: Atoms,
    vib_kwargs: dict = None,
    scratch_dir: str = SETTINGS.SCRATCH_DIR,
    gzip: bool = SETTINGS.GZIP_FILES,
    copy_files: list[str] = None,
) -> Atoms:
    """
    Run an ASE-based vibration analysis in a scratch directory and copy the results
    back to the original directory. This can be useful if file I/O is slow in
    the working directory, so long as file transfer speeds are reasonable.

    This is a wrapper around the vibrations module in ASE. Note: This function does
    not modify the atoms object in-place.

    Parameters
    ----------
    atoms
        The Atoms object to run the calculation on.
    vib_kwargs
        Dictionary of kwargs for the vibration analysis.
    scratch_dir
        Path where a tmpdir should be made for running the calculation. If None,
        the current working directory will be used.
    gzip
        Whether to gzip the output files.
    copy_files
        Filenames to copy from source to scratch directory.

    Returns
    -------
    .Vibrations
        The updated Vibrations module
    """

    if atoms.calc is None:
        raise ValueError("Atoms object must have attached calculator.")

    atoms = copy_atoms(atoms)
    cwd = os.getcwd()
    scratch_dir = scratch_dir or cwd
    symlink = os.path.join(cwd, "tmp_dir")
    vib_kwargs = vib_kwargs or {}

    if not os.path.exists(scratch_dir):
        os.makedirs(scratch_dir)

    tmpdir = mkdtemp(prefix="quacc-tmp", dir=scratch_dir)

    if os.name != "nt":
        if os.path.islink(symlink):
            os.unlink(symlink)
        os.symlink(tmpdir, symlink)

    # Copy files to scratch and decompress them if needed
    if copy_files:
        copy_decompress(copy_files, tmpdir)

    # Run calculation
    os.chdir(tmpdir)
    vib = Vibrations(atoms, **vib_kwargs)
    vib.run()
    vib.summary()
    os.chdir(cwd)

    # Gzip files in tmpdir
    if gzip:
        gzip_dir(tmpdir)

    # Copy files back to run_dir
    copy_r(tmpdir, cwd)

    # Remove symlink
    if os.path.islink(symlink):
        os.remove(symlink)

    os.chdir(cwd)

    return vib


def ideal_gas_thermo(
    atoms: Atoms,
    vib_list: list[float | complex],
    temperature: float = 298.15,
    pressure: float = 1.0,
    energy: float = 0.0,
    spin_multiplicity: float = None,
) -> dict:
    """
    Calculate thermodynamic properties for a molecule from a given vibrational analysis.

    Parameters
    ----------
    atoms
        The Atoms object to use.
    vib_list
        The list of vibrations to use, typically obtained from Vibrations.get_frequencies().
    temperature
        Temperature in Kelvins.
    pressure
        Pressure in bar.
    energy
        Potential energy in eV. If 0 eV, then the thermochemical correction is computed.
    spin_multiplicity
        The spin multiplicity. If None, this will be determined automatically from the
        attached magnetic moments.

    Returns
    -------
    dict
        {
            "atoms": .Atoms,
            ...,
            "results":
            {
                "frequencies": list of frequencies in cm^-1,
                "true_frequencies": list of true vibrational frequencies in cm^-1,
                "n_imag": number of imaginary modes within true_frequencies in cm^-1,
                "geometry": the geometry of the molecule,
                "pointgroup": the point group of the molecule,
                "energy": potential energy in eV,
                "enthalpy": enthalpy in eV,
                "entropy": entropy in eV/K,
                "gibbs_energy": free energy in eV
            }
        }
    """
    for i, f in enumerate(vib_list):
        if not isinstance(f, complex) and f < 0:
            vib_list[i] = complex(0 - f * 1j)

    # Get the spin from the Atoms object
    if spin_multiplicity:
        spin = (spin_multiplicity - 1) / 2
    elif (
        getattr(atoms, "calc", None) is not None
        and getattr(atoms.calc, "results", None) is not None
    ):
        spin = round(atoms.calc.results.get("magmom", 0)) / 2
    elif atoms.has("initial_magmoms"):
        spin = round(np.sum(atoms.get_initial_magnetic_moments())) / 2
    else:
        spin = 0

    # Get symmetry for later use
    natoms = len(atoms)
    pmg_obj = AseAtomsAdaptor.get_molecule(atoms)
    pga = PointGroupAnalyzer(pmg_obj)

    pointgroup = None if len(atoms) == 1 else pga.get_pointgroup().sch_symbol

    # Get the geometry and true frequencies that should
    # be used for thermo calculations
    if natoms == 1:
        geometry = "monatomic"
        true_freqs = []
    elif natoms == 2 or (natoms > 2 and pointgroup == "D*h"):
        geometry = "linear"
        true_freqs = vib_list[-(3 * natoms - 5) :]
    else:
        geometry = "nonlinear"
        true_freqs = vib_list[-(3 * natoms - 6) :]

    # Automatically get rotational symmetry number
    if geometry == "monatomic":
        symmetry_number = 1
    else:
        symmetry_number = pga.get_rotational_symmetry_number()

    # Fetch the real vibrational energies
    real_vib_freqs = [f for f in true_freqs if np.isreal(f)]
    n_imag = len(true_freqs) - len(real_vib_freqs)
    real_vib_energies = [f * invcm for f in real_vib_freqs]

    # Calculate ideal gas thermo
    igt = IdealGasThermo(
        real_vib_energies,
        geometry,
        potentialenergy=energy,
        atoms=atoms,
        symmetrynumber=symmetry_number,
        spin=spin,
    )
    if len(igt.vib_energies) != len(real_vib_energies):
        raise ValueError(
            "The number of real vibrational modes and those used by ASE do not match. Something is very wrong..."
        )

    # Use negative sign convention for imag modes
    vib_list = [np.abs(f) if np.isreal(f) else -np.abs(f) for f in vib_list]
    true_freqs = [np.abs(f) if np.isreal(f) else -np.abs(f) for f in true_freqs]

    return {
        **atoms_to_metadata(atoms),
        "results": {
            "frequencies": vib_list,  # full list of computed frequencies
            "true_frequencies": true_freqs,  # list of *relevant* frequencies based on the geometry
            "n_imag": n_imag,  # number of imag modes within true_frequencies
            "geometry": geometry,
            "pointgroup": pointgroup,
            "energy": energy,
            "enthalpy": igt.get_enthalpy(temperature, verbose=False),
            "entropy": igt.get_entropy(temperature, pressure * 10**5, verbose=False),
            "gibbs_energy": igt.get_gibbs_energy(
                temperature, pressure * 10**5, verbose=False
            ),
        },
    }


def _check_logfile(logfile: str, check_str: str) -> bool:
    """
    Check if a logfile has a given string (case-insensitive).

    Parameters
    ----------
    logfile : str
        Path to the logfile.
    check_str : str
        String to check for.

    Returns
    -------
    bool
        True if the string is found in the logfile, False otherwise.
    """
    zlog = zpath(logfile)
    with zopen(zlog, "r") as f:
        for line in f:
            if not isinstance(line, str):
                line = line.decode("utf-8")
            if check_str.lower() in line.lower():
                return True
    return False


def _convert_auto_kpts(
    atoms: Atoms,
    auto_kpts: None
    | dict[str, float]
    | dict[str, list[tuple[float, float]]]
    | dict[str, list[tuple[float, float, float]]],
    force_gamma: bool = True,
) -> tuple[list[tuple[int, int, int]], None | bool, None | bool]:
    """
    Shortcuts for pymatgen k-point generation schemes.
    Options include: line_density (for band structures),
    reciprocal_density (by volume), grid_density (by number of atoms),
    max_mixed_density (max of reciprocal_density volume or atoms), and length_density (good for slabs).
    These are formatted as {"line_density": float}, {"reciprocal_density": float},
    {"grid_density": float}, {"max_mixed_density": [float, float]}, and
    {"length_density": [float, float, float]}.

    Parameters
    ----------
    .Atoms
        ASE Atoms object
    auto_kpts
        Dictionary describing the automatic k-point scheme
    force_gamma
        Whether a gamma-centered mesh should be returned

    Returns
    -------
    List[int, int, int]
        List of k-points for use with the ASE Vasp calculator
    Optional[bool]
        The gamma command for use with the ASE Vasp calculator
    Optional[bool]
        The reciprocal command for use with the ASE Vasp calculator
    """
    struct = AseAtomsAdaptor.get_structure(atoms)

    if auto_kpts.get("line_density", None):
        # TODO: Support methods other than latimer-munro
        kpath = HighSymmKpath(struct, path_type="latimer_munro")
        kpts, _ = kpath.get_kpoints(
            line_density=auto_kpts["line_density"], coords_are_cartesian=True
        )
        kpts = np.stack(kpts)
        reciprocal = True
        gamma = None

    else:
        reciprocal = None
        if auto_kpts.get("max_mixed_density", None):
            if len(auto_kpts["max_mixed_density"]) != 2:
                raise ValueError("Must specify two values for max_mixed_density.")

            if auto_kpts["max_mixed_density"][0] > auto_kpts["max_mixed_density"][1]:
                warnings.warn(
                    "Warning: It is not usual that kppvol > kppa. Please make sure you have chosen the right k-point densities.",
                )
            pmg_kpts1 = Kpoints.automatic_density_by_vol(
                struct, auto_kpts["max_mixed_density"][0], force_gamma=force_gamma
            )
            pmg_kpts2 = Kpoints.automatic_density(
                struct, auto_kpts["max_mixed_density"][1], force_gamma=force_gamma
            )
            if np.product(pmg_kpts1.kpts[0]) >= np.product(pmg_kpts2.kpts[0]):
                pmg_kpts = pmg_kpts1
            else:
                pmg_kpts = pmg_kpts2
        elif auto_kpts.get("reciprocal_density", None):
            pmg_kpts = Kpoints.automatic_density_by_vol(
                struct, auto_kpts["reciprocal_density"], force_gamma=force_gamma
            )
        elif auto_kpts.get("grid_density", None):
            pmg_kpts = Kpoints.automatic_density(
                struct, auto_kpts["grid_density"], force_gamma=force_gamma
            )
        elif auto_kpts.get("length_density", None):
            if len(auto_kpts["length_density"]) != 3:
                raise ValueError("Must specify three values for length_density.")
            pmg_kpts = Kpoints.automatic_density_by_lengths(
                struct, auto_kpts["length_density"], force_gamma=force_gamma
            )
        else:
            raise ValueError(f"Unsupported k-point generation scheme: {auto_kpts}.")

        kpts = pmg_kpts.kpts[0]
        gamma = pmg_kpts.style.name.lower() == "gamma"

    return kpts, gamma, reciprocal
