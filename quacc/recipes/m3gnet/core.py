"""
Core recipes for the MatGL code
"""
from __future__ import annotations

from copy import deepcopy

import covalent as ct
from ase.atoms import Atoms
from ase.constraints import ExpCellFilter
from monty.dev import requires

from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.util.calc import run_ase_opt, run_calc

try:
    from matgl.ext.ase import M3GNetCalculator as M3GNet
except ImportError:
    M3GNet = None
