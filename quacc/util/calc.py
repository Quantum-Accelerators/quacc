"""
Utility functions for running ASE calculators
"""
from __future__ import annotations

import os
from shutil import copy
from tempfile import TemporaryDirectory
from typing import Any, Dict, List

from ase.atoms import Atoms
from ase.io import read
from ase.optimize import FIRE
from ase.optimize.optimize import Optimizer
from monty.os.path import zpath
from monty.shutil import copy_r, decompress_file, gzip_dir

from quacc import SETTINGS
from quacc.util.atoms import copy_atoms


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

    with TemporaryDirectory(dir=scratch_dir) as tmpdir:

        if os.name != "nt":
            os.symlink(tmpdir, symlink)

        # Copy files to scratch and decompress them if needed
        if copy_files:
            for f in copy_files:
                zpath = zpath(f)
                zfile = os.path.basename(zpath)
                if os.path.exists(zpath):
                    copy(zpath, os.path.join(tmpdir, zfile))
                    decompress_file(os.path.join(tmpdir, zfile))

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

    with TemporaryDirectory(dir=scratch_dir) as tmpdir:

        if os.name != "nt":
            os.symlink(tmpdir, symlink)

        # Copy files to scratch and decompress them if needed
        if copy_files:
            for f in copy_files:
                zpath = zpath(f)
                zfile = os.path.basename(zpath)
                if os.path.exists(zpath):
                    copy(zpath, os.path.join(tmpdir, zfile))
                    decompress_file(os.path.join(tmpdir, zfile))

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
