from __future__ import annotations

from importlib.util import find_spec
from pathlib import Path
from shutil import copy

import pytest
from ase import units
from ase.calculators.lj import LennardJones
from ase.io import read
from ase.optimize import FIRE
from pymatgen.io.qchem.inputs import QCInput

from quacc import _internally_set_settings
from quacc.atoms.core import check_charge_and_spin
from quacc.calculators.qchem import QChem
from quacc.recipes.qchem.core import freq_job, relax_job, static_job
from quacc.recipes.qchem.ts import irc_job, quasi_irc_job, ts_job

has_sella = bool(find_spec("sella"))


FILE_DIR = Path(__file__).parent
QCHEM_DIR = FILE_DIR / "qchem_examples"


def setup_module():
    _internally_set_settings({"CHECK_CONVERGENCE": False})


def teardown_module():
    _internally_set_settings(reset=True)


@pytest.fixture()
def test_atoms():
    return read(FILE_DIR / "xyz" / "test.xyz")


@pytest.fixture()
def test_qirc_atoms():
    return read(FILE_DIR / "xyz" / "ts_test.xyz")


@pytest.fixture()
def os_atoms():
    return read(FILE_DIR / "xyz" / "OS_test.xyz")


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


def mock_execute1(self, **kwargs):
    copy(QCHEM_DIR / "mol.qout.basic", Path(self.directory, "mol.qout"))
    copy(QCHEM_DIR / "131.0.basic", Path(self.directory, "131.0"))
    copy(QCHEM_DIR / "53.0.basic", Path(self.directory, "53.0"))
    copy(QCHEM_DIR / "custodian.json", Path(self.directory, "custodian.json"))


def mock_execute2(self, **kwargs):
    copy(QCHEM_DIR / "mol.qout.intermediate", Path(self.directory, "mol.qout"))
    copy(QCHEM_DIR / "131.0.intermediate", Path(self.directory, "131.0"))
    copy(QCHEM_DIR / "53.0.intermediate", Path(self.directory, "53.0"))
    copy(QCHEM_DIR / "custodian.json", Path(self.directory, "custodian.json"))


def mock_execute3(self, **kwargs):
    copy(QCHEM_DIR / "mol.qout.alternate", Path(self.directory, "mol.qout"))
    copy(QCHEM_DIR / "131.0.alternate", Path(self.directory, "131.0"))
    copy(QCHEM_DIR / "53.0.alternate", Path(self.directory, "53.0"))
    copy(QCHEM_DIR / "custodian.json", Path(self.directory, "custodian.json"))


def mock_execute4(self, **kwargs):
    qcin = QCInput.from_file(str(Path(self.directory, "mol.qin")))
    mol = qcin.molecule
    atoms = mol.to_ase_atoms()
    atoms.calc = LennardJones()
    atoms.get_potential_energy()
    self.results = atoms.calc.results


def mock_execute5(self, **kwargs):
    copy(QCHEM_DIR / "mol.qout.freq", Path(self.directory, "mol.qout"))
    copy(QCHEM_DIR / "132.0.freq", Path(self.directory, "132.0"))
    copy(QCHEM_DIR / "53.0.freq", Path(self.directory, "53.0"))
    copy(QCHEM_DIR / "custodian.json", Path(self.directory, "custodian.json"))


def mock_read(self, **kwargs):
    if self.results is None:
        raise RuntimeError("Results should not be None here.")


def test_static_job_v1(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(QChem, "execute", mock_execute1)
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

    qcin = QCInput.from_file(str(Path(output["dir_name"], "mol.qin.gz")))
    ref_qcin = QCInput.from_file(str(QCHEM_DIR / "mol.qin.basic"))
    qcinput_nearly_equal(qcin, ref_qcin)
    assert output["results"]["taskdoc"]


def test_static_job_v2(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)

    monkeypatch.setattr(QChem, "execute", mock_execute2)
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

    qcin = QCInput.from_file(str(Path(output["dir_name"], "mol.qin.gz")))
    ref_qcin = QCInput.from_file(str(QCHEM_DIR / "mol.qin.intermediate"))
    qcinput_nearly_equal(qcin, ref_qcin)
    assert output["results"]["taskdoc"]


def test_static_job_v3(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)

    monkeypatch.setattr(QChem, "execute", mock_execute3)
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

    qcin = QCInput.from_file(str(Path(output["dir_name"], "mol.qin.gz")))
    ref_qcin = QCInput.from_file(str(QCHEM_DIR / "mol.qin.alternate"))
    qcinput_nearly_equal(qcin, ref_qcin)
    assert output["results"]["taskdoc"]


def test_static_job_v4(monkeypatch, tmp_path, os_atoms):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(QChem, "read_results", mock_read)
    monkeypatch.setattr(QChem, "execute", mock_execute4)
    charge, spin_multiplicity = check_charge_and_spin(os_atoms)
    assert static_job(os_atoms, charge=charge, spin_multiplicity=spin_multiplicity)


def test_static_job_v5(tmp_path, monkeypatch, test_atoms):
    monkeypatch.chdir(tmp_path)

    with pytest.raises(
        ValueError,
        match="Only one of PCM, ISOSVP, SMD, and CMIRSmay be used for solvation",
    ):
        static_job(
            test_atoms,
            charge=0,
            spin_multiplicity=1,
            qchem_dict_set_params={"pcm_dielectric": "3.0", "smd_solvent": "water"},
        )


@pytest.mark.skipif(has_sella is False, reason="Does not have Sella")
def test_relax_job_v1(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)

    monkeypatch.setattr(QChem, "execute", mock_execute1)
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

    qcin = QCInput.from_file(str(Path(output["dir_name"], "mol.qin.gz")))
    ref_qcin = QCInput.from_file(str(QCHEM_DIR / "mol.qin.basic.sella_opt_iter1"))
    qcinput_nearly_equal(qcin, ref_qcin)
    assert len(output["results"]["taskdoc"]["input"]) > 1


@pytest.mark.skipif(has_sella is False, reason="Does not have Sella")
def test_relax_job_v2(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(QChem, "execute", mock_execute2)
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

    qcin = QCInput.from_file(str(Path(output["dir_name"], "mol.qin.gz")))
    ref_qcin = QCInput.from_file(
        str(QCHEM_DIR / "mol.qin.intermediate.sella_opt_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)


@pytest.mark.skipif(has_sella is False, reason="Does not have Sella")
def test_relax_job_v3(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(QChem, "execute", mock_execute3)
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


@pytest.mark.skipif(has_sella is False, reason="Does not have Sella")
def test_relax_job_v4(tmp_path, monkeypatch, test_atoms):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(
        ValueError,
        match="Only one of PCM, ISOSVP, SMD, and CMIRSmay be used for solvation",
    ):
        relax_job(
            test_atoms,
            charge=0,
            spin_multiplicity=1,
            qchem_dict_set_params={"pcm_dielectric": "3.0", "smd_solvent": "water"},
        )


def test_freq_job_v1(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(QChem, "execute", mock_execute5)
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
    assert output["results"]["taskdoc"]["output"]["enthalpy"] is not None


@pytest.mark.skipif(has_sella is False, reason="Does not have Sella")
def test_ts_job_v1(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)

    monkeypatch.setattr(QChem, "execute", mock_execute1)
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

    qcin = QCInput.from_file(str(Path(output["dir_name"], "mol.qin.gz")))
    ref_qcin = QCInput.from_file(str(QCHEM_DIR / "mol.qin.basic.sella_TSopt_iter1"))
    qcinput_nearly_equal(qcin, ref_qcin)


@pytest.mark.skipif(has_sella is False, reason="Does not have Sella")
def test_ts_job_v2(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(QChem, "execute", mock_execute2)
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

    qcin = QCInput.from_file(str(Path(output["dir_name"], "mol.qin.gz")))
    ref_qcin = QCInput.from_file(
        str(QCHEM_DIR / "mol.qin.intermediate.sella_TSopt_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)


@pytest.mark.skipif(has_sella is False, reason="Does not have Sella")
def test_ts_job_v3(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(QChem, "execute", mock_execute3)
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


@pytest.mark.skipif(has_sella is False, reason="Does not have Sella")
def test_ts_job_v4(tmp_path, monkeypatch, test_atoms):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(
        ValueError,
        match="Only one of PCM, ISOSVP, SMD, and CMIRSmay be used for solvation",
    ):
        ts_job(
            test_atoms,
            charge=0,
            spin_multiplicity=1,
            qchem_dict_set_params={"pcm_dielectric": "3.0", "smd_solvent": "water"},
        )

    with pytest.raises(
        ValueError, match="Only Sella should be used for TS optimization"
    ):
        ts_job(
            test_atoms,
            charge=0,
            spin_multiplicity=1,
            qchem_dict_set_params={"pcm_dielectric": "3.0", "smd_solvent": "water"},
            opt_params={"optimizer": FIRE},
        )


@pytest.mark.skipif(has_sella is False, reason="Does not have Sella")
def test_irc_job_v1(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)

    monkeypatch.setattr(QChem, "read_results", mock_read)
    monkeypatch.setattr(QChem, "execute", mock_execute4)

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

    qcin = QCInput.from_file(str(Path(output["dir_name"], "mol.qin.gz")))
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

    qcin = QCInput.from_file(str(Path(output["dir_name"], "mol.qin.gz")))
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


@pytest.mark.skipif(has_sella is False, reason="Does not have Sella")
def test_irc_job_v2(tmp_path, monkeypatch, test_atoms):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(
        ValueError, match='direction must be one of "forward" or "reverse"!'
    ):
        irc_job(test_atoms, charge=0, spin_multiplicity=1, direction="straight")

    with pytest.raises(
        ValueError,
        match="Only one of PCM, ISOSVP, SMD, and CMIRSmay be used for solvation",
    ):
        irc_job(
            test_atoms,
            charge=0,
            spin_multiplicity=1,
            direction="forward",
            qchem_dict_set_params={"pcm_dielectric": "3.0", "smd_solvent": "water"},
        )

    with pytest.raises(
        ValueError, match="Only Sella's IRC should be used for IRC optimization"
    ):
        irc_job(
            test_atoms,
            charge=0,
            spin_multiplicity=1,
            direction="forward",
            qchem_dict_set_params={"pcm_dielectric": "3.0", "smd_solvent": "water"},
            opt_params={"optimizer": FIRE},
        )


@pytest.mark.skipif(has_sella is False, reason="Does not have Sella")
def test_quasi_irc_job(monkeypatch, tmp_path, test_qirc_atoms):
    monkeypatch.chdir(tmp_path)

    monkeypatch.setattr(QChem, "read_results", mock_read)
    monkeypatch.setattr(QChem, "execute", mock_execute4)

    # Transition mode for this transition-state
    mode = [
        [-0.164, 0.289, 0.027],
        [0.112, -0.02, -0.004],
        [0.012, -0.072, -0.042],
        [-0.087, 0.039, -0.038],
        [-0.017, 0.013, 0.001],
        [0.028, -0.186, 0.028],
        [0.751, -0.378, 0.186],
        [0.042, 0.034, 0.025],
        [-0.007, -0.001, -0.009],
        [-0.056, -0.179, -0.076],
        [0.036, 0.035, 0.027],
        [0.043, 0.037, 0.023],
        [0.036, 0.032, 0.021],
        [-0.003, -0.032, 0.011],
        [-0.006, -0.009, -0.118],
        [0.014, -0.034, 0.094],
    ]

    charge, spin_multiplicity = check_charge_and_spin(test_qirc_atoms)

    output = quasi_irc_job(
        test_qirc_atoms,
        mode,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        direction="forward",
        method="wb97mv",
        opt_params={"max_steps": 5},
        basis="def2-svpd",
    )

    assert output["atoms"] != test_qirc_atoms
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H8 O4"
    assert output["nelectrons"] == 64
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1

    qcin = QCInput.from_file(str(Path(output["dir_name"], "mol.qin.gz")))
    ref_qcin = QCInput.from_file(str(QCHEM_DIR / "mol.qin.qirc_forward"))
    qcinput_nearly_equal(qcin, ref_qcin)

    output = quasi_irc_job(
        test_qirc_atoms,
        mode,
        charge=-1,
        spin_multiplicity=2,
        perturb_magnitude=1.0,
        direction="reverse",
        basis="def2-tzvpd",
        opt_params={"max_steps": 6},
        rem={"scf_algorithm": "gdm"},
    )

    assert output["atoms"] != test_qirc_atoms
    assert output["charge"] == -1
    assert output["spin_multiplicity"] == 2
    assert output["formula_alphabetical"] == "C4 H8 O4"
    assert output["nelectrons"] == 65
    assert output["parameters"]["charge"] == -1
    assert output["parameters"]["spin_multiplicity"] == 2

    qcin = QCInput.from_file(str(Path(output["dir_name"], "mol.qin.gz")))
    ref_qcin = QCInput.from_file(str(QCHEM_DIR / "mol.qin.qirc_reverse"))
    qcinput_nearly_equal(qcin, ref_qcin)
