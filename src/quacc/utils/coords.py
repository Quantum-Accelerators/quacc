"""Utilities for working with atomic coordinates."""

from __future__ import annotations
import copy

import numpy as np

from ase import Atoms


def perturb(mol: Atoms, vector: List[List[float]], scale: float):
    mol_copy = copy.deepcopy(mol)
    mode_copy = copy.deepcopy(vector)

    orig_pos = mol_copy.get_positions()

    if not isinstance(mode_copy, np.ndarray):
        mode_copy = np.asarray(mode_copy)

    pos = orig_pos + mode_copy * scale
    mol_copy.set_positions(pos)

    return mol_copy