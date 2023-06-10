"""
Core recipes for the Sella code
"""
from __future__ import annotations

from copy import deepcopy

import covalent as ct
import numpy as np
from ase.atoms import Atoms
from ase.io import read
from monty.dev import requires

try:
    from sella import Sella
except ImportError:
    Sella = None
try:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet
except ImportError:
    NewtonNet = None

# TODO: Don't hardcode fmax or steps
# TODO: Add docstrings
# TODO: Add all path-related variables to the global Quacc settings and remove from kwargs.


@ct.electron
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
@requires(Sella, "Sella must be installed. Try pip install quacc[sella]")
def optimize_and_analyze(atoms: Atoms, ml_path: str, config_path: str) -> dict:
    mlcalculator = NewtonNet(model_path=ml_path, settings_path=config_path)
    atoms.calc = mlcalculator
    opt = Sella(atoms, internal=True, logfile=f"sella.log", trajectory=f"sella.traj")
    opt.run(fmax=0.01, steps=5)

    traj = read(f"sella.traj", index=":")
    mlcalculator.calculate(traj)
    H = mlcalculator.results["hessian"]
    n_atoms = np.shape(H)[0]
    A = np.reshape(H, (n_atoms * 3, n_atoms * 3))
    eigvals, eigvecs = np.linalg.eig(A)
    return_vars = {
        "eigvals": np.array2string(eigvals),
        "eigvecs": np.array2string(eigvecs),
    }

    summary = {
        "name": name,
        "input": {"atoms": deepcopy(atoms)},
        "output": {"return_vars": return_vars},
    }
    return summary
