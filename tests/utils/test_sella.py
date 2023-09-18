import pytest
from ase.build import bulk, molecule
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones

from quacc.utils.calc import run_ase_opt

try:
    from sella import Sella
except ImportError:
    Sella = None


@pytest.mark.skipif(
    not Sella,
    reason="Sella must be installed.",
)
def test_sella(tmpdir):
    tmpdir.chdir()

    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()
    dyn = run_ase_opt(atoms, optimizer=Sella, optimizer_kwargs={"restart": None})
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None
    assert dyn.user_internal is False

    atoms = molecule("H2O")
    atoms.calc = LennardJones()
    dyn = run_ase_opt(atoms, optimizer=Sella, optimizer_kwargs={"restart": None})
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None
    assert dyn.user_internal is True


@pytest.mark.skipif(
    not Sella,
    reason="Sella must be installed.",
)
def test_TRICs(tmpdir):
    tmpdir.chdir()
    atoms = molecule("O2")
    atoms.calc = LennardJones()
    dyn = run_ase_opt(
        atoms,
        optimizer=Sella,
        optimizer_kwargs={"use_TRICs": True},
    )
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None
