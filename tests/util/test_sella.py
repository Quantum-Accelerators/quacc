import os

import pytest
from ase.build import bulk, molecule
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones

from quacc.util.calc import run_ase_opt

try:
    import sella
except ImportError:
    sella = None


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_sella(tmpdir):
    tmpdir.chdir()
    from sella.optimize import Sella

    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()
    dyn = run_ase_opt(
        atoms,
        optimizer=Sella,
        scratch_dir="test_calc",
        gzip=False,
        optimizer_kwargs={"restart": None},
    )
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None
    assert dyn.user_internal is False

    atoms = molecule("H2O")
    atoms.calc = LennardJones()
    dyn = run_ase_opt(
        atoms,
        optimizer=Sella,
        scratch_dir="test_calc2",
        gzip=False,
        optimizer_kwargs={"restart": None},
    )
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None
    assert dyn.user_internal is True
