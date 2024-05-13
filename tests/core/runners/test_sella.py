from __future__ import annotations

import pytest

sella = pytest.importorskip("sella")
from ase.build import bulk, molecule
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones
from sella import Sella

from quacc.runners.ase import run_opt


def test_sella():
    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()
    dyn = run_opt(atoms, optimizer=Sella, optimizer_kwargs={"restart": None})
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None
    assert dyn.user_internal is False

    atoms = molecule("H2O")
    atoms.calc = LennardJones()
    dyn = run_opt(atoms, optimizer=Sella, optimizer_kwargs={"restart": None})
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None
    assert dyn.user_internal is True


def test_dof():
    atoms = molecule("C2H6")
    atoms.calc = LennardJones()
    dyn = run_opt(atoms, optimizer=Sella, fmax=1.0)
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None
    assert dyn.internal.ndof == 24
    assert dyn.internal.nbonds == 0
    assert dyn.internal.nangles == 0
    assert dyn.internal.ndihedrals == 0

    atoms = bulk("Cu")
    atoms.calc = EMT()
    dyn = run_opt(atoms, optimizer=Sella, fmax=1.0)
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None
    assert dyn.internal is None
