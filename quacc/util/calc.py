"""
Utility functions for running ASE calculators
"""
from __future__ import annotations

import os
import warnings
from tempfile import mkdtemp
from typing import Any, Dict, List, Tuple

import numpy as np
from ase.atoms import Atoms
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
from ase.vibrations import Vibrations
from ase.units import invcm
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
from quacc.util.basics import copy_decompress


def run_calc(
    atoms: Atoms,
    geom_file: str = None,
    scratch_dir: str = SETTINGS.SCRATCH_DIR,
    gzip: bool = SETTINGS.GZIP_FILES,
    copy_files: List[str] = None,
) -> Atoms:
    """
    Run a calculation in a scratch directory and copy the results back to the
    original directory. This can be useful if file I/O is slow in the working
    directory, so long as file transfer speeds are reasonable.

    This is a wrapper around atoms.get_potential_energy(). Note: This
    function does not modify the atoms object in-place.

    Parameters
    ----------
    atoms : .Atoms
        The Atoms object to run the calculation on.
    geom_file : str
        The filename of the log file that contains the output geometry, used
        to update the atoms object's positions and cell after a job. It is better
        to specify this rather than relying on ASE's atoms.get_potential_energy()
        function to update the positions, as this varies between codes.
    scratch_dir : str
        Path where a tmpdir should be made for running the calculation. If None,
        the current working directory will be used.
    gzip : bool
        Whether to gzip the output files.
    copy_files : List[str]
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
        atoms.positions = atoms_new.positions
        atoms.cell = atoms_new.cell

    return atoms


def run_ase_opt(
    atoms: Atoms,
    fmax: float = 0.01,
    max_steps: int = 100,
    optimizer: str = "FIRE",
    opt_kwargs: Dict[str, Any] = None,
    scratch_dir: str = SETTINGS.SCRATCH_DIR,
    gzip: bool = SETTINGS.GZIP_FILES,
    copy_files: List[str] = None,
) -> trajectory:
    """
    Run an ASE-based optimization in a scratch directory and copy the results
    back to the original directory. This can be useful if file I/O is slow in
    the working directory, so long as file transfer speeds are reasonable.

    This is a wrapper around the optimizers in ASE. Note: This function does
    not modify the atoms object in-place.

    Parameters
    ----------
    atoms : .Atoms
        The Atoms object to run the calculation on.
    fmax : float
        Tolerance for the force convergence (in eV/A).
    max_steps : int
        Maximum number of steps to take.
    optimizer : str
        Name of optimizer class to use.
    opt_kwargs : dict
        Dictionary of kwargs for the optimizer.
    scratch_dir : str
        Path where a tmpdir should be made for running the calculation. If None,
        the current working directory will be used.
    gzip : bool
        Whether to gzip the output files.
    copy_files : List[str]
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
    vib_kwargs: Dict[str, Any] = None,
    scratch_dir: str = SETTINGS.SCRATCH_DIR,
    gzip: bool = SETTINGS.GZIP_FILES,
    copy_files: List[str] = None,
) -> Atoms:
    """
    Run an ASE-based vibration analysis in a scratch directory and copy the results
    back to the original directory. This can be useful if file I/O is slow in
    the working directory, so long as file transfer speeds are reasonable.

    This is a wrapper around the vibrations module in ASE. Note: This function does
    not modify the atoms object in-place.

    Parameters
    ----------
    atoms : .Atoms
        The Atoms object to run the calculation on.
    vib_kwargs : dict
        Dictionary of kwargs for the vibration analysis.
    scratch_dir : str
        Path where a tmpdir should be made for running the calculation. If None,
        the current working directory will be used.
    gzip : bool
        Whether to gzip the output files.
    copy_files : List[str]
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
    vibrations: Vibrations,
    atoms: Atoms = None,
    temperature: float = 298.15,
    pressure: float = 1.0,
    energy: float = 0.0,
    spin_multiplicity: float = None,
) -> Dict[str, Any]:
    """
    Calculate thermodynamic properties for a molecule from a given vibrational analysis.

    Parameters
    ----------
    vibrations : .Vibrations
        The Vibrations module to use.
    atoms : .Atoms
        The Atoms object to use. If None, the Atoms object will be taken from
        the Vibrations module.
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
                "n_imag": number of imaginary modes,
                "energy": potential energy in eV,
                "enthalpy": enthalpy in eV,
                "entropy": entropy in eV/K,
                "gibbs_energy": free energy in eV
            }
        }
    """
    # Pull atoms from vibrations object if needed
    atoms = atoms or vibrations.atoms

    # Get the spin from the Atoms object
    if spin_multiplicity:
        spin = (spin_multiplicity - 1) / 2
    else:
        if (
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

    if len(atoms) == 1:
        pointgroup = None
    else:
        pointgroup = pga.get_pointgroup().sch_symbol

    # Get the geometry and true frequencies that should
    # be used for thermo calculations
    all_freqs = vibrations.get_frequencies()
    if natoms == 1:
        geometry = "monatomic"
        true_freqs = []
    elif natoms == 2 or (natoms > 2 and pointgroup == "D*h"):
        geometry = "linear"
        true_freqs = all_freqs[-(3 * natoms - 5) :]
    else:
        geometry = "nonlinear"
        true_freqs = all_freqs[-(3 * natoms - 6) :]

    # Automatically get rotational symmetry number
    if geometry == "monatomic":
        symmetry_number = 1
    else:
        symmetry_number = pga.get_rotational_symmetry_number()

    # Fetch the real vibrational energies
    real_vib_freqs = [f for f in true_freqs if np.isreal(f)]
    n_imag = len(true_freqs) - len(real_vib_freqs)
    real_vib_energies = [f * invcm for f in true_freqs]

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
    all_freqs = [np.abs(f) if np.isreal(f) else -np.abs(f) for f in all_freqs]
    true_freqs = [np.abs(f) if np.isreal(f) else -np.abs(f) for f in true_freqs]

    # Count number of relevant imaginary frequencies

    thermo_summary = {
        **atoms_to_metadata(atoms),
        "results": {
            "frequencies": all_freqs,  # full list of computed frequencies
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

    return thermo_summary


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
    | Dict[str, float]
    | Dict[str, List[Tuple[float, float]]]
    | Dict[str, List[Tuple[float, float, float]]],
    force_gamma: bool = True,
) -> Tuple[List[Tuple[int, int, int]], None | bool, None | bool]:
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
