import os
from pathlib import Path
from shutil import copy, rmtree

from ase.build import molecule
from pymatgen.core import Molecule
from pymatgen.io.ase import AseAtomsAdaptor

from quacc.recipes.qchem.core import relax_job, static_job

FILE_DIR = Path(__file__).resolve().parent
QCHEM_DIR = os.path.join(FILE_DIR, "qchem_run")


def setup_module():
    for f in os.listdir(QCHEM_DIR):
        copy(os.path.join(QCHEM_DIR, f), os.path.join(os.getcwd(), f))


def teardown_module():
    for f in os.listdir(QCHEM_DIR):
        if os.path.exists(os.path.join(os.getcwd(), f)):
            os.remove(os.path.join(os.getcwd(), f))
    for f in os.listdir(os.getcwd()):
        if "quacc-tmp" in f or f == "tmp_dir":
            if os.path.islink(f):
                os.unlink(f)
            else:
                rmtree(f)


def test_static_Job():
    mol = Molecule.from_file(os.path.join(FILE_DIR, "test.xyz"))
    atoms = AseAtomsAdaptor.get_atoms(mol)

    output = static_job(atoms)
    assert output["atoms"] == atoms
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == 'C4 H4 O6'
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == None
    assert output["parameters"]["spin_multiplicity"] == None
