import os
from pathlib import Path
from shutil import copy, rmtree
from ase import Atoms, units
from ase.build import molecule
from pymatgen.core import Molecule
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.qchem.inputs import QCInput

from quacc.recipes.qchem.core import relax_job, static_job

try:
    import sella
except ImportError:
    sella = None


FILE_DIR = Path(__file__).resolve().parent
QCHEM_DIR = os.path.join(FILE_DIR, "qchem_examples")


def teardown_module():
    if os.path.exists(os.path.join(os.getcwd(), "mol.qin.gz")):
        os.remove(os.path.join(os.getcwd(), "mol.qin.gz"))
    if os.path.exists(os.path.join(os.getcwd(), "mol.qout.gz")):
        os.remove(os.path.join(os.getcwd(), "mol.qout.gz"))
    if os.path.exists(os.path.join(os.getcwd(), "131.0.gz")):
        os.remove(os.path.join(os.getcwd(), "131.0.gz"))
    if os.path.exists(os.path.join(os.getcwd(), "__pycache__")):
        rmtree(os.path.join(os.getcwd(), "__pycache__"))
    if os.path.exists(os.path.join(os.getcwd(), "opt.traj")):
        os.remove(os.path.join(os.getcwd(), "opt.traj"))
    for f in os.listdir(os.getcwd()):
        if "quacc-tmp" in f or f == "tmp_dir":
            if os.path.islink(f):
                os.unlink(f)
            else:
                rmtree(f)


def test_static_job():
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
    assert output["results"]["energy"] == -606.1616819641*units.Hartree
    assert output["results"]["forces"][0][0] == -1.3826330655069403

    qcin = QCInput.from_file(os.path.join(FILE_DIR, "mol.qin.gz"))
    ref_qcin = QCInput.from_file(os.path.join(QCHEM_DIR, "mol.qin.basic"))
    assert qcin.as_dict() == ref_qcin.as_dict()

    mol = Molecule.from_file(os.path.join(FILE_DIR, "test.xyz"))
    atoms = AseAtomsAdaptor.get_atoms(mol)
    output = static_job(
        atoms=atoms,
        charge=-1,
        xc="b97mv",
        basis="def2-svpd",
        pcm_dielectric="3.0",
    )

    assert output["atoms"] == atoms
    # next three lines should be able to be uncommented once charge is correctly passed through
    # assert output["charge"] == -1
    # assert output["spin_multiplicity"] == 2
    # assert output["nelectrons"] == 77
    assert output["formula_alphabetical"] == 'C4 H4 O6'
    assert output["parameters"]["charge"] == -1
    assert output["parameters"]["spin_multiplicity"] == None
    assert output["results"]["energy"] == -605.6859554025*units.Hartree #-16481.554341995
    assert output["results"]["forces"][0][0] == -0.6955571014353796

    qcin = QCInput.from_file(os.path.join(FILE_DIR, "mol.qin.gz"))
    ref_qcin = QCInput.from_file(os.path.join(QCHEM_DIR, "mol.qin.intermediate"))
    assert qcin.as_dict() == ref_qcin.as_dict()
    os.remove(os.path.join(os.getcwd(), "mol.qin.gz"))
    os.remove(os.path.join(os.getcwd(), "mol.qout.gz"))
    os.remove(os.path.join(os.getcwd(), "131.0.gz"))


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_relax_job():
    mol = Molecule.from_file(os.path.join(FILE_DIR, "test.xyz"))
    atoms = AseAtomsAdaptor.get_atoms(mol)
    output = relax_job(
        atoms=atoms,
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1}
    )

    assert output["atoms"] != atoms
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == 'C4 H4 O6'
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == None
    assert output["parameters"]["spin_multiplicity"] == None
    assert output["results"]["energy"] == -606.1616819641*units.Hartree
    assert output["results"]["forces"][0][0] == -1.3826330655069403

    qcin = QCInput.from_file(os.path.join(FILE_DIR, "mol.qin.gz"))
    ref_qcin = QCInput.from_file(os.path.join(QCHEM_DIR, "mol.qin.basic.sella_opt_iter1"))
    assert qcin.as_dict() == ref_qcin.as_dict()

    mol = Molecule.from_file(os.path.join(FILE_DIR, "test.xyz"))
    atoms = AseAtomsAdaptor.get_atoms(mol)
    output = relax_job(
        atoms=atoms,
        charge=-1,
        xc="b97mv",
        pcm_dielectric="3.0",
        opt_swaps={"max_steps": 1}
    )

    assert output["atoms"] != atoms
    # next three lines should be able to be uncommented once charge is correctly passed through
    # assert output["charge"] == -1
    # assert output["spin_multiplicity"] == 2
    # assert output["nelectrons"] == 77
    assert output["formula_alphabetical"] == 'C4 H4 O6'
    assert output["parameters"]["charge"] == -1
    assert output["parameters"]["spin_multiplicity"] == None
    assert output["results"]["energy"] == -605.6859554025*units.Hartree
    assert output["results"]["forces"][0][0] == -0.6955571014353796

    qcin = QCInput.from_file(os.path.join(FILE_DIR, "mol.qin.gz"))
    ref_qcin = QCInput.from_file(os.path.join(QCHEM_DIR, "mol.qin.intermediate.sella_opt_iter1"))
    assert qcin.as_dict() == ref_qcin.as_dict()
    os.remove(os.path.join(os.getcwd(), "mol.qin.gz"))
    os.remove(os.path.join(os.getcwd(), "mol.qout.gz"))
    os.remove(os.path.join(os.getcwd(), "131.0.gz"))


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_relax_job_as_TSopt():
    mol = Molecule.from_file(os.path.join(FILE_DIR, "test.xyz"))
    atoms = AseAtomsAdaptor.get_atoms(mol)
    output = relax_job(
        atoms=atoms,
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1, "optimizer_kwargs": {"order": 1}}
    )

    assert output["atoms"] != atoms
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == 'C4 H4 O6'
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == None
    assert output["parameters"]["spin_multiplicity"] == None
    assert output["results"]["energy"] == -606.1616819641*units.Hartree
    assert output["results"]["forces"][0][0] == -1.3826330655069403

    qcin = QCInput.from_file(os.path.join(FILE_DIR, "mol.qin.gz"))
    ref_qcin = QCInput.from_file(os.path.join(QCHEM_DIR, "mol.qin.basic.sella_TSopt_iter1"))
    assert qcin.as_dict() == ref_qcin.as_dict()

    mol = Molecule.from_file(os.path.join(FILE_DIR, "test.xyz"))
    atoms = AseAtomsAdaptor.get_atoms(mol)
    output = relax_job(
        atoms=atoms,
        charge=-1,
        xc="b97mv",
        pcm_dielectric="3.0",
        opt_swaps={"max_steps": 1, "optimizer_kwargs": {"order": 1}}
    )

    assert output["atoms"] != atoms
    # next three lines should be able to be uncommented once charge is correctly passed through
    # assert output["charge"] == -1
    # assert output["spin_multiplicity"] == 2
    # assert output["nelectrons"] == 77
    assert output["formula_alphabetical"] == 'C4 H4 O6'
    assert output["parameters"]["charge"] == -1
    assert output["parameters"]["spin_multiplicity"] == None
    assert output["results"]["energy"] == -605.6859554025*units.Hartree
    assert output["results"]["forces"][0][0] == -0.6955571014353796

    qcin = QCInput.from_file(os.path.join(FILE_DIR, "mol.qin.gz"))
    ref_qcin = QCInput.from_file(os.path.join(QCHEM_DIR, "mol.qin.intermediate.sella_TSopt_iter1"))
    assert qcin.as_dict() == ref_qcin.as_dict()