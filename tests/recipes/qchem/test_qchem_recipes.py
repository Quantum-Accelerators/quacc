import os
from pathlib import Path
from shutil import copy, rmtree

import pytest
from ase import units
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.lj import LennardJones
from ase.io import read
from ase.optimize import FIRE
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.qchem.inputs import QCInput

from quacc.calculators.qchem import QChem
from quacc.recipes.qchem.core import (
    irc_job,
    quasi_irc_job,
    relax_job,
    static_job,
    ts_job,
)

try:
    import sella
except ImportError:
    sella = None


FILE_DIR = Path(__file__).resolve().parent
QCHEM_DIR = os.path.join(FILE_DIR, "qchem_examples")
TEST_ATOMS = read(os.path.join(FILE_DIR, "test.xyz"))
OS_ATOMS = read(os.path.join(FILE_DIR, "OS_test.xyz"))


def qcinput_nearly_equal(qcinput1, qcinput2):
    qcin1 = qcinput1.as_dict()
    qcin2 = qcinput2.as_dict()
    for key in qcin1:
        if key == "molecule":
            for molkey in qcin1[key]:
                if molkey == "sites":
                    for ii, site in enumerate(qcin1[key][molkey]):
                        for sitekey in site:
                            if sitekey == "xyz":
                                for jj, val in enumerate(site[sitekey]):
                                    assert val == pytest.approx(
                                        qcin2[key][molkey][ii][sitekey][jj]
                                    )
                            else:
                                assert (
                                    qcin1[key][molkey][ii][sitekey]
                                    == qcin2[key][molkey][ii][sitekey]
                                )

                else:
                    assert qcin1[key][molkey] == qcin2[key][molkey]

        else:
            assert qcin1[key] == qcin2[key]


def mock_execute1(_self, **kwargs):
    copy(os.path.join(QCHEM_DIR, "mol.qout.basic"), "mol.qout")
    copy(os.path.join(QCHEM_DIR, "131.0.basic"), "131.0")
    copy(os.path.join(QCHEM_DIR, "53.0.basic"), "53.0")


def mock_execute2(_self, **kwargs):
    copy(os.path.join(QCHEM_DIR, "mol.qout.intermediate"), "mol.qout")
    copy(os.path.join(QCHEM_DIR, "131.0.intermediate"), "131.0")
    copy(os.path.join(QCHEM_DIR, "53.0.intermediate"), "53.0")


def mock_execute3(_self, **kwargs):
    copy(os.path.join(QCHEM_DIR, "mol.qout.alternate"), "mol.qout")
    copy(os.path.join(QCHEM_DIR, "131.0.alternate"), "131.0")
    copy(os.path.join(QCHEM_DIR, "53.0.alternate"), "53.0")


def mock_execute4(self, **kwargs):
    qcin = QCInput.from_file("mol.qin")
    mol = qcin.molecule
    atoms = AseAtomsAdaptor.get_atoms(mol)
    atoms.calc = LennardJones()
    atoms.get_potential_energy()
    self.results = atoms.calc.results


def mock_read(self, **kwargs):
    if self.results is None:
        raise RuntimeError("Results should not be None here!")
    pass


def teardown_module():
    for f in os.listdir("."):
        if ".log" in f or ".traj" in f or ".gz" in f:
            os.remove(f)
    for f in os.listdir(os.getcwd()):
        if "quacc-tmp" in f or f == "tmp_dir":
            if os.path.islink(f):
                os.unlink(f)
            else:
                rmtree(f)


def test_static_job(monkeypatch):
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute1)
    output = static_job(TEST_ATOMS)
    assert output["atoms"] == TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] is None
    assert output["parameters"]["spin_multiplicity"] is None
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826330655069403)

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(os.path.join(QCHEM_DIR, "mol.qin.basic"))
    qcinput_nearly_equal(qcin, ref_qcin)

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute2)
    output = static_job(
        atoms=TEST_ATOMS,
        charge=-1,
        method="b97mv",
        basis="def2-svpd",
        pcm_dielectric="3.0",
    )

    assert output["atoms"] == TEST_ATOMS
    assert output["charge"] == -1
    assert output["spin_multiplicity"] == 2
    assert output["nelectrons"] == 77
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["parameters"]["charge"] == -1
    assert output["parameters"]["spin_multiplicity"] is None
    assert output["results"]["energy"] == pytest.approx(-605.6859554025 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-0.6955571014353796)

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(os.path.join(QCHEM_DIR, "mol.qin.intermediate"))
    qcinput_nearly_equal(qcin, ref_qcin)

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute3)
    overwrite_inputs = {"rem": {"mem_total": "170000"}}
    output = static_job(
        TEST_ATOMS, scf_algorithm="gdm", overwrite_inputs=overwrite_inputs
    )
    assert output["atoms"] == TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] is None
    assert output["parameters"]["spin_multiplicity"] is None
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826311086011256)

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(os.path.join(QCHEM_DIR, "mol.qin.alternate"))
    qcinput_nearly_equal(qcin, ref_qcin)

    with pytest.raises(ValueError):
        output = static_job(atoms=TEST_ATOMS, pcm_dielectric="3.0", smd_solvent="water")

    monkeypatch.setattr(QChem, "read_results", mock_read)
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute4)
    output = static_job(OS_ATOMS)


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_relax_job(monkeypatch):
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute1)
    output = relax_job(
        atoms=TEST_ATOMS,
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
        check_convergence=False,
    )

    assert output["atoms"] != TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] is None
    assert output["parameters"]["spin_multiplicity"] is None
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826330655069403)

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(
        os.path.join(QCHEM_DIR, "mol.qin.basic.sella_opt_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute2)
    output = relax_job(
        atoms=TEST_ATOMS,
        charge=-1,
        method="b97mv",
        pcm_dielectric="3.0",
        opt_swaps={"max_steps": 1},
        check_convergence=False,
    )

    assert output["atoms"] != TEST_ATOMS
    assert output["charge"] == -1
    assert output["spin_multiplicity"] == 2
    assert output["nelectrons"] == 77
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["parameters"]["charge"] == -1
    assert output["parameters"]["spin_multiplicity"] is None
    assert output["results"]["energy"] == pytest.approx(-605.6859554025 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-0.6955571014353796)

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(
        os.path.join(QCHEM_DIR, "mol.qin.intermediate.sella_opt_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute3)
    overwrite_inputs = {"rem": {"mem_total": "170000"}}
    output = relax_job(
        atoms=TEST_ATOMS,
        scf_algorithm="gdm",
        overwrite_inputs=overwrite_inputs,
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
        check_convergence=False,
    )

    assert output["atoms"] != TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] is None
    assert output["parameters"]["spin_multiplicity"] is None
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826311086011256)

    with pytest.raises(ValueError):
        output = relax_job(atoms=TEST_ATOMS, pcm_dielectric="3.0", smd_solvent="water")


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_ts_job(monkeypatch):
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute1)
    output = ts_job(
        atoms=TEST_ATOMS,
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
        check_convergence=False,
    )

    assert output["atoms"] != TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] is None
    assert output["parameters"]["spin_multiplicity"] is None
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826330655069403)

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(
        os.path.join(QCHEM_DIR, "mol.qin.basic.sella_TSopt_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute2)
    output = ts_job(
        atoms=TEST_ATOMS,
        charge=-1,
        method="b97mv",
        pcm_dielectric="3.0",
        opt_swaps={"max_steps": 1},
        check_convergence=False,
    )

    assert output["atoms"] != TEST_ATOMS
    assert output["charge"] == -1
    assert output["spin_multiplicity"] == 2
    assert output["nelectrons"] == 77
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["parameters"]["charge"] == -1
    assert output["parameters"]["spin_multiplicity"] is None
    assert output["results"]["energy"] == pytest.approx(-605.6859554025 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-0.6955571014353796)

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(
        os.path.join(QCHEM_DIR, "mol.qin.intermediate.sella_TSopt_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute3)
    overwrite_inputs = {"rem": {"mem_total": "170000"}}
    output = ts_job(
        atoms=TEST_ATOMS,
        scf_algorithm="gdm",
        overwrite_inputs=overwrite_inputs,
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
        check_convergence=False,
    )

    assert output["atoms"] != TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] is None
    assert output["parameters"]["spin_multiplicity"] is None
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826311086011256)

    with pytest.raises(ValueError):
        output = ts_job(atoms=TEST_ATOMS, pcm_dielectric="3.0", smd_solvent="water")

    with pytest.raises(RuntimeError):
        output = ts_job(
            atoms=TEST_ATOMS,
            pcm_dielectric="3.0",
            smd_solvent="water",
            opt_swaps={"optimizer": FIRE},
        )


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_irc_job(monkeypatch):
    monkeypatch.setattr(QChem, "read_results", mock_read)
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute4)

    output = irc_job(
        atoms=TEST_ATOMS,
        direction="forward",
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
        check_convergence=False,
    )

    assert output["atoms"] != TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] is None
    assert output["parameters"]["spin_multiplicity"] is None

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(
        os.path.join(QCHEM_DIR, "mol.qin.basic.sella_IRC_forward_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)

    output = irc_job(
        atoms=TEST_ATOMS,
        direction="reverse",
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
        check_convergence=False,
    )

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(
        os.path.join(QCHEM_DIR, "mol.qin.basic.sella_IRC_reverse_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)

    overwrite_inputs = {"rem": {"mem_total": "170000"}}
    output = irc_job(
        atoms=TEST_ATOMS,
        direction="reverse",
        scf_algorithm="gdm",
        overwrite_inputs=overwrite_inputs,
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
        check_convergence=False,
    )

    assert output["atoms"] != TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] is None
    assert output["parameters"]["spin_multiplicity"] is None

    with pytest.raises(ValueError):
        output = irc_job(atoms=TEST_ATOMS, direction="straight")

    with pytest.raises(ValueError):
        output = irc_job(
            atoms=TEST_ATOMS,
            direction="forward",
            pcm_dielectric="3.0",
            smd_solvent="water",
        )

    with pytest.raises(ValueError):
        output = irc_job(
            atoms=TEST_ATOMS,
            direction="forward",
            pcm_dielectric="3.0",
            smd_solvent="water",
            opt_swaps={"optimizer": FIRE},
        )


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_quasi_irc_job(monkeypatch):
    monkeypatch.setattr(QChem, "read_results", mock_read)
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute4)

    common_kwargs = {"basis": "def2-tzvpd"}
    relax_opt_swaps = {"max_steps": 5}

    output = quasi_irc_job(
        atoms=TEST_ATOMS,
        direction="forward",
        common_kwargs=common_kwargs,
        relax_opt_swaps=relax_opt_swaps,
        check_convergence=False,
    )

    assert output["atoms"] != TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] is None
    assert output["parameters"]["spin_multiplicity"] is None

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(
        os.path.join(QCHEM_DIR, "mol.qin.basic.quasi_irc_forward")
    )
    qcinput_nearly_equal(qcin, ref_qcin)

    common_kwargs = {"charge": -1, "basis": "def2-svpd", "scf_algorithm": "gdm"}
    irc_opt_swaps = {"max_steps": 6}
    relax_opt_swaps = {"max_steps": 6}

    output = quasi_irc_job(
        atoms=TEST_ATOMS,
        direction="reverse",
        common_kwargs=common_kwargs,
        irc_opt_swaps=irc_opt_swaps,
        relax_opt_swaps=relax_opt_swaps,
        check_convergence=False,
    )

    assert output["atoms"] != TEST_ATOMS
    assert output["charge"] == -1
    assert output["spin_multiplicity"] == 2
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 77
    assert output["parameters"]["charge"] == -1
    assert output["parameters"]["spin_multiplicity"] is None

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(os.path.join(QCHEM_DIR, "mol.qin.quasi_irc_reverse"))
    qcinput_nearly_equal(qcin, ref_qcin)
