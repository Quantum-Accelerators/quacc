"""
Core recipes for the Sella code
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import Any, Dict

import covalent as ct
import numpy as np
from ase.atoms import Atoms
from ase.io import read
from jobflow import Maker, job
from monty.dev import requires

try:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet
except ImportError:
    NewtonNet = None
from pymatgen.io.ase import AseAtomsAdaptor
from sella import Sella


@dataclass
@requires(NewtonNet is not None, "NewtonNet must be installed")
class OptimizeAndAnalyzeJob(Maker):
    name: str = "Optimize-And-Analyze"

    @job
    def make(self, atoms: Atoms, ml_path: str, config_path: str) -> Dict[str, Any]:
        ase_atoms_adaptor = AseAtomsAdaptor()
        atoms = ase_atoms_adaptor.get_atoms(atoms)

        mlcalculator = NewtonNet(model_path=ml_path, settings_path=config_path)
        atoms.set_calculator(mlcalculator)
        opt = Sella(
            atoms, internal=True, logfile=f"sella.log", trajectory=f"sella.traj"
        )
        opt.run(fmax=0.01, steps=5)

        traj = read(f"sella.traj")
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
