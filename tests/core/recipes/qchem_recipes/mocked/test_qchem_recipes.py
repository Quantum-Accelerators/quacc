from pathlib import Path
from shutil import copy

import pytest
from ase import units
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.lj import LennardJones
from ase.io import read
from ase.optimize import FIRE
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.qchem.inputs import QCInput

from quacc import SETTINGS
from quacc.atoms.core import check_charge_and_spin
from quacc.calculators.qchem import QChem
from quacc.recipes.qchem.core import freq_job, relax_job, static_job
from quacc.recipes.qchem.ts import irc_job, quasi_irc_job, ts_job

try:
    import sella
except ImportError:
    sella = None

FILE_DIR = Path(__file__).parent
QCHEM_DIR = FILE_DIR / "qchem_examples"

DEFAULT_SETTINGS = SETTINGS.model_copy()


@pytest.fixture()
def test_atoms():
    return read(FILE_DIR / "xyz" / "test.xyz")


@pytest.fixture()
def os_atoms():
    return read(FILE_DIR / "xyz" / "OS_test.xyz")


def setup_module():
    SETTINGS.CHECK_CONVERGENCE = False


def teardown_module():
    SETTINGS.CHECK_CONVERGENCE = DEFAULT_SETTINGS.CHECK_CONVERGENCE


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
    copy(QCHEM_DIR / "mol.qout.basic", "mol.qout")
    copy(QCHEM_DIR / "131.0.basic", "131.0")
    copy(QCHEM_DIR / "53.0.basic", "53.0")
    copy(QCHEM_DIR / "custodian.json", "custodian.json")


def mock_execute2(_self, **kwargs):
    copy(QCHEM_DIR / "mol.qout.intermediate", "mol.qout")
    copy(QCHEM_DIR / "131.0.intermediate", "131.0")
    copy(QCHEM_DIR / "53.0.intermediate", "53.0")
    copy(QCHEM_DIR / "custodian.json", "custodian.json")


def mock_execute3(_self, **kwargs):
    copy(QCHEM_DIR / "mol.qout.alternate", "mol.qout")
    copy(QCHEM_DIR / "131.0.alternate", "131.0")
    copy(QCHEM_DIR / "53.0.alternate", "53.0")
    copy(QCHEM_DIR / "custodian.json", "custodian.json")


def mock_execute4(self, **kwargs):
    qcin = QCInput.from_file("mol.qin")
    mol = qcin.molecule
    atoms = AseAtomsAdaptor.get_atoms(mol)
    atoms.calc = LennardJones()
    atoms.get_potential_energy()
    self.results = atoms.calc.results


def mock_execute5(_self, **kwargs):
    copy(QCHEM_DIR / "mol.qout.freq", "mol.qout")
    copy(QCHEM_DIR / "132.0.freq", "132.0")
    copy(QCHEM_DIR / "53.0.freq", "53.0")
    copy(QCHEM_DIR / "custodian.json", "custodian.json")


def mock_read(self, **kwargs):
    if self.results is None:
        raise RuntimeError("Results should not be None here.")


def test_static_job_v1(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute1)
    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = static_job(test_atoms, charge=charge, spin_multiplicity=spin_multiplicity)
    assert output["atoms"] == test_atoms
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826330655069403)
    assert output["results"]["custodian"][0]["job"]["max_cores"] == 40

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(str(QCHEM_DIR / "mol.qin.basic"))
    qcinput_nearly_equal(qcin, ref_qcin)
    qcinput_nearly_equal(ref_qcin, QCInput.from_dict(output["results"]["qc_input"]))


def test_static_job_v2(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute2)
    charge, spin_multiplicity = check_charge_and_spin(test_atoms, charge=-1)
    output = static_job(
        test_atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        method="b97mv",
        basis="def2-svpd",
        qchem_dict_set_params={"pcm_dielectric": "3.0"},
    )

    assert output["atoms"] == test_atoms
    assert output["charge"] == -1
    assert output["spin_multiplicity"] == 2
    assert output["nelectrons"] == 77
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["parameters"]["charge"] == -1
    assert output["parameters"]["spin_multiplicity"] == 2
    assert output["results"]["energy"] == pytest.approx(-605.6859554025 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-0.6955571014353796)

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(str(QCHEM_DIR / "mol.qin.intermediate"))
    qcinput_nearly_equal(qcin, ref_qcin)
    qcinput_nearly_equal(ref_qcin, QCInput.from_dict(output["results"]["qc_input"]))


def test_static_job_v3(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute3)
    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = static_job(
        test_atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        rem={"mem_total": 170000, "scf_algorithm": "gdm"},
    )
    assert output["atoms"] == test_atoms
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826311086011256)

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(str(QCHEM_DIR / "mol.qin.alternate"))
    qcinput_nearly_equal(qcin, ref_qcin)
    qcinput_nearly_equal(ref_qcin, QCInput.from_dict(output["results"]["qc_input"]))


def test_static_job_v4(monkeypatch, tmp_path, os_atoms):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(QChem, "read_results", mock_read)
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute4)
    charge, spin_multiplicity = check_charge_and_spin(os_atoms)
    assert static_job(os_atoms, charge=charge, spin_multiplicity=spin_multiplicity)


def test_static_job_v5(tmp_path, monkeypatch, test_atoms):
    monkeypatch.chdir(tmp_path)

    with pytest.raises(ValueError):
        static_job(
            test_atoms,
            charge=0,
            spin_multiplicity=1,
            qchem_dict_set_params={"pcm_dielectric": "3.0", "smd_solvent": "water"},
        )


@pytest.mark.skipif(sella is None, reason="Does not have Sella")
def test_relax_job_v1(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute1)
    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = relax_job(
        test_atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        basis="def2-tzvpd",
        opt_params={"max_steps": 1},
    )

    assert output["atoms"] != test_atoms
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826330655069403)

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(str(QCHEM_DIR / "mol.qin.basic.sella_opt_iter1"))
    qcinput_nearly_equal(qcin, ref_qcin)
    assert len(output["results"]["qc_input"]) > 1


@pytest.mark.skipif(sella is None, reason="Does not have Sella")
def test_relax_job_v2(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute2)
    charge, spin_multiplicity = check_charge_and_spin(test_atoms, charge=-1)
    output = relax_job(
        test_atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        method="b97mv",
        qchem_dict_set_params={"pcm_dielectric": "3.0"},
        opt_params={"max_steps": 1},
    )

    assert output["atoms"] != test_atoms
    assert output["charge"] == -1
    assert output["spin_multiplicity"] == 2
    assert output["nelectrons"] == 77
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["parameters"]["charge"] == -1
    assert output["parameters"]["spin_multiplicity"] == 2
    assert output["results"]["energy"] == pytest.approx(-605.6859554025 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-0.6955571014353796)

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(
        str(QCHEM_DIR / "mol.qin.intermediate.sella_opt_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)


@pytest.mark.skipif(sella is None, reason="Does not have Sella")
def test_relax_job_v3(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute3)
    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = relax_job(
        test_atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        rem={"scf_algorithm": "gdm", "mem_total": 170000},
        basis="def2-tzvpd",
        opt_params={"max_steps": 1},
    )

    assert output["atoms"] != test_atoms
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826311086011256)


@pytest.mark.skipif(sella is None, reason="Does not have Sella")
def test_relax_job_v4(tmp_path, monkeypatch, test_atoms):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError):
        relax_job(
            test_atoms,
            charge=0,
            spin_multiplicity=1,
            qchem_dict_set_params={"pcm_dielectric": "3.0", "smd_solvent": "water"},
        )


def test_freq_job_v1(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute5)
    charge, spin_multiplicity = check_charge_and_spin(test_atoms, charge=-1)
    output = freq_job(
        test_atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        rem={"scf_algorithm": "diis"},
        method="b97mv",
        basis="def2-svpd",
    )

    assert output["atoms"] == test_atoms
    assert output["charge"] == -1
    assert output["spin_multiplicity"] == 2
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 77
    assert output["parameters"]["charge"] == -1
    assert output["parameters"]["spin_multiplicity"] == 2
    assert output["results"]["energy"] == pytest.approx(-605.6859554019 * units.Hartree)
    assert output["results"].get("hessian") is not None
    assert output["results"]["enthalpy"] == pytest.approx(
        output["results"]["qc_output"]["total_enthalpy"] * (units.kcal / units.mol)
    )
    assert (
        output["results"]["entropy"]
        == output["results"]["qc_output"]["total_entropy"]
        * 0.001
        * units.kcal
        / units.mol
    )


@pytest.mark.skipif(sella is None, reason="Does not have Sella")
def test_ts_job_v1(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute1)
    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = ts_job(
        test_atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        basis="def2-tzvpd",
        opt_params={"max_steps": 1},
    )

    assert output["atoms"] != test_atoms
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826330655069403)

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(str(QCHEM_DIR / "mol.qin.basic.sella_TSopt_iter1"))
    qcinput_nearly_equal(qcin, ref_qcin)


@pytest.mark.skipif(sella is None, reason="Does not have Sella")
def test_ts_job_v2(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute2)
    charge, spin_multiplicity = check_charge_and_spin(test_atoms, charge=-1)
    output = ts_job(
        test_atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        method="b97mv",
        qchem_dict_set_params={"pcm_dielectric": "3.0"},
        opt_params={"max_steps": 1},
    )

    assert output["atoms"] != test_atoms
    assert output["charge"] == -1
    assert output["spin_multiplicity"] == 2
    assert output["nelectrons"] == 77
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["parameters"]["charge"] == -1
    assert output["parameters"]["spin_multiplicity"] == 2
    assert output["results"]["energy"] == pytest.approx(-605.6859554025 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-0.6955571014353796)

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(
        str(QCHEM_DIR / "mol.qin.intermediate.sella_TSopt_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)


@pytest.mark.skipif(sella is None, reason="Does not have Sella")
def test_ts_job_v3(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute3)
    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = ts_job(
        test_atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        rem={"scf_algorithm": "gdm", "mem_total": 170000},
        basis="def2-tzvpd",
        opt_params={"max_steps": 1},
    )

    assert output["atoms"] != test_atoms
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826311086011256)


@pytest.mark.skipif(sella is None, reason="Does not have Sella")
def test_ts_job_v4(tmp_path, monkeypatch, test_atoms):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError):
        ts_job(
            test_atoms,
            charge=0,
            spin_multiplicity=1,
            qchem_dict_set_params={"pcm_dielectric": "3.0", "smd_solvent": "water"},
        )

    with pytest.raises(ValueError):
        ts_job(
            test_atoms,
            charge=0,
            spin_multiplicity=1,
            qchem_dict_set_params={"pcm_dielectric": "3.0", "smd_solvent": "water"},
            opt_params={"optimizer": FIRE},
        )


@pytest.mark.skipif(sella is None, reason="Does not have Sella")
def test_irc_job_v1(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)

    monkeypatch.setattr(QChem, "read_results", mock_read)
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute4)

    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = irc_job(
        test_atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        direction="forward",
        basis="def2-tzvpd",
        opt_params={"max_steps": 1},
    )

    assert output["atoms"] != test_atoms
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(
        str(QCHEM_DIR / "mol.qin.basic.sella_IRC_forward_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)

    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = irc_job(
        test_atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        direction="reverse",
        basis="def2-tzvpd",
        opt_params={"max_steps": 1},
    )

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(
        str(QCHEM_DIR / "mol.qin.basic.sella_IRC_reverse_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)

    output = irc_job(
        test_atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        direction="reverse",
        rem={"scf_algorithm": "gdm", "mem_total": 170000},
        basis="def2-tzvpd",
        opt_params={"max_steps": 1},
    )

    assert output["atoms"] != test_atoms
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1


@pytest.mark.skipif(sella is None, reason="Does not have Sella")
def test_irc_job_v2(tmp_path, monkeypatch, test_atoms):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError):
        irc_job(test_atoms, charge=0, spin_multiplicity=1, direction="straight")

    with pytest.raises(ValueError):
        irc_job(
            test_atoms,
            charge=0,
            spin_multiplicity=1,
            direction="forward",
            qchem_dict_set_params={"pcm_dielectric": "3.0", "smd_solvent": "water"},
        )

    with pytest.raises(ValueError):
        irc_job(
            test_atoms,
            charge=0,
            spin_multiplicity=1,
            direction="forward",
            qchem_dict_set_params={"pcm_dielectric": "3.0", "smd_solvent": "water"},
            opt_params={"optimizer": FIRE},
        )


@pytest.mark.skipif(sella is None, reason="Does not have Sella")
def test_quasi_irc_job(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)

    monkeypatch.setattr(QChem, "read_results", mock_read)
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute4)

    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = quasi_irc_job(
        test_atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        direction="forward",
        basis="def2-tzvpd",
        relax_job_kwargs={"opt_params": {"max_steps": 5}},
    )

    assert output["atoms"] != test_atoms
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(str(QCHEM_DIR / "mol.qin.basic.quasi_irc_forward"))
    qcinput_nearly_equal(qcin, ref_qcin)

    output = quasi_irc_job(
        test_atoms,
        charge=-1,
        spin_multiplicity=2,
        direction="reverse",
        basis="def2-svpd",
        irc_job_kwargs={
            "rem": {"scf_algorithm": "gdm"},
            "opt_params": {"max_steps": 6},
        },
        relax_job_kwargs={
            "rem": {"scf_algorithm": "gdm"},
            "opt_params": {"max_steps": 6},
        },
    )

    assert output["atoms"] != test_atoms
    assert output["charge"] == -1
    assert output["spin_multiplicity"] == 2
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 77
    assert output["parameters"]["charge"] == -1
    assert output["parameters"]["spin_multiplicity"] == 2

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(str(QCHEM_DIR / "mol.qin.quasi_irc_reverse"))
    qcinput_nearly_equal(qcin, ref_qcin)
