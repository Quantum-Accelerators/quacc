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
from ase.io import read
from ase.optimize import FIRE
from ase.optimize.optimize import Optimizer
from ase.thermochemistry import IdealGasThermo
from ase.vibrations import Vibrations
from monty.io import zopen
from monty.os.path import zpath
from monty.shutil import copy_r, gzip_dir
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath

from quacc import SETTINGS
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

    tmpdir = mkdtemp(dir=scratch_dir)

    if os.name != "nt":
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
    optimizer: Optimizer = FIRE,
    opt_kwargs: Dict[str, Any] = None,
    scratch_dir: str = SETTINGS.SCRATCH_DIR,
    gzip: bool = SETTINGS.GZIP_FILES,
    copy_files: List[str] = None,
) -> Atoms:
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
    optimizer : .Optimizer
        .Optimizer class to use for the relaxation.
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
    .Atoms
        The updated .Atoms object,
    """

    if atoms.calc is None:
        raise ValueError("Atoms object must have attached calculator.")

    atoms = copy_atoms(atoms)
    cwd = os.getcwd()
    scratch_dir = scratch_dir or cwd
    symlink = os.path.join(cwd, "tmp_dir")
    opt_kwargs = opt_kwargs or {}

    tmpdir = mkdtemp(dir=scratch_dir)

    if os.name != "nt":
        os.symlink(tmpdir, symlink)

    # Copy files to scratch and decompress them if needed
    if copy_files:
        copy_decompress(copy_files, tmpdir)

    # Run calculation
    os.chdir(tmpdir)
    dyn = optimizer(atoms, **opt_kwargs)
    dyn.run(fmax=fmax)
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

    return atoms


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

    tmpdir = mkdtemp(dir=scratch_dir)

    if os.name != "nt":
        os.symlink(tmpdir, symlink)

    # Copy files to scratch and decompress them if needed
    if copy_files:
        copy_decompress(copy_files, tmpdir)

    # Run calculation
    os.chdir(tmpdir)
    vib = Vibrations(atoms, **vib_kwargs)
    vib.run()
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


def calculate_thermo(
    vibrations: Vibrations,
    atoms: Atoms = None,
    temperature: float = 298.15,
    pressure: float = 1.0,
    energy: float = 0.0,
    geometry: str = None,
    symmetry_number: int = 1,
    spin_multiplicity: float = None,
) -> Dict[str, Any]:
    """
    Calculate thermodynamic properties from a given vibrational analysis.

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
    geometry
        Monatomic, linear, or nonlinear. Will try to determine automatically if None.
    symmetry_number
        Rotational symmetry number.
    spin_multiplicity
        The spin multiplicity

    Returns
    -------
    dict
        {"frequencies": list of frequencies in cm^-1,
        "enthalpy": enthalpy in eV,
        "entropy": entropy in eV/K,
        "free_energy": free energy in eV}
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
    if geometry is None or symmetry_number is None:
        if atoms.pbc.any():
            pmg_obj = AseAtomsAdaptor.get_structure(atoms)
        else:
            pmg_obj = AseAtomsAdaptor.get_molecule(atoms)
        pga = PointGroupAnalyzer(pmg_obj)
        pointgroup = pga.get_pointgroup()

    # Get the geometry
    if geometry is None:
        if len(atoms) == 1:
            geometry = "monatomic"
        elif len(atoms) == 2 or (len(atoms) > 2 and pointgroup == "D*h"):
            geometry = "linear"
        else:
            geometry = "nonlinear"

    # TODO: Automatically get rotational symmetry number if None

    # Calculate ideal gas thermo
    igt = IdealGasThermo(
        vibrations.get_energies(),
        geometry,
        potentialenergy=energy,
        atoms=atoms,
        symmetrynumber=symmetry_number,
        spin=spin,
    )
    freqs = vibrations.get_frequencies()

    # Use negataive sign convention for imag modes
    clean_freqs = []
    for f in freqs:
        if np.iscomplex(f):
            clean_freqs.append(-np.abs(f))
        else:
            clean_freqs.append(np.abs(f))

    thermo_summary = {
        "frequencies": clean_freqs,
        "enthalpy": igt.get_enthalpy(temperature, verbose=False),
        "entropy": igt.get_entropy(temperature, pressure / 10**5, verbose=False),
        "free_energy": igt.get_gibbs_energy(
            temperature, pressure / 10**5, verbose=False
        ),
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
            if check_str.lower() in line.decode("utf-8").lower():
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
