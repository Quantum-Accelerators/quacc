import os
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
from quacc.calculators.qchem import QChem
from quacc.recipes.qchem.core import freq_job, internal_relax_job, relax_job, static_job
from quacc.recipes.qchem.ts import irc_job, quasi_irc_job, ts_job
from quacc.utils import check_charge_and_spin

try:
    import sella
except ImportError:
    sella = None


FILE_DIR = Path(__file__).resolve().parent
QCHEM_DIR = os.path.join(FILE_DIR, "qchem_examples")
TEST_ATOMS = read(os.path.join(FILE_DIR, "test.xyz"))
OS_ATOMS = read(os.path.join(FILE_DIR, "OS_test.xyz"))

DEFAULT_SETTINGS = SETTINGS.copy()


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


def mock_execute5(_self, **kwargs):
    copy(os.path.join(QCHEM_DIR, "mol.qout.freq"), "mol.qout")
    copy(os.path.join(QCHEM_DIR, "132.0.freq"), "132.0")
    copy(os.path.join(QCHEM_DIR, "53.0.freq"), "53.0")


def mock_read(self, **kwargs):
    if self.results is None:
        raise RuntimeError("Results should not be None here.")


def test_static_job_v1(monkeypatch, tmpdir):
    tmpdir.chdir()

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute1)
    charge, spin_multiplicity = check_charge_and_spin(TEST_ATOMS)
    output = static_job(TEST_ATOMS, charge, spin_multiplicity)
    assert output["atoms"] == TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826330655069403)

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(os.path.join(QCHEM_DIR, "mol.qin.basic"))
    qcinput_nearly_equal(qcin, ref_qcin)


def test_static_job_v2(monkeypatch, tmpdir):
    tmpdir.chdir()

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute2)
    charge, spin_multiplicity = check_charge_and_spin(TEST_ATOMS, charge=-1)
    output = static_job(
        TEST_ATOMS,
        charge,
        spin_multiplicity,
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
    assert output["parameters"]["spin_multiplicity"] == 2
    assert output["results"]["energy"] == pytest.approx(-605.6859554025 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-0.6955571014353796)

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(os.path.join(QCHEM_DIR, "mol.qin.intermediate"))
    qcinput_nearly_equal(qcin, ref_qcin)


def test_static_job_v3(monkeypatch, tmpdir):
    tmpdir.chdir()

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute3)
    overwrite_inputs = {"rem": {"mem_total": "170000"}}
    charge, spin_multiplicity = check_charge_and_spin(TEST_ATOMS)
    output = static_job(
        TEST_ATOMS,
        charge,
        spin_multiplicity,
        scf_algorithm="gdm",
        overwrite_inputs=overwrite_inputs,
    )
    assert output["atoms"] == TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826311086011256)

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(os.path.join(QCHEM_DIR, "mol.qin.alternate"))
    qcinput_nearly_equal(qcin, ref_qcin)


def test_static_job_v4(monkeypatch, tmpdir):
    tmpdir.chdir()
    monkeypatch.setattr(QChem, "read_results", mock_read)
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute4)
    charge, spin_multiplicity = check_charge_and_spin(OS_ATOMS)
    assert static_job(OS_ATOMS, charge, spin_multiplicity)


def test_static_job_v5(tmpdir):
    tmpdir.chdir()

    with pytest.raises(ValueError):
        static_job(TEST_ATOMS, 0, 1, pcm_dielectric="3.0", smd_solvent="water")


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_relax_job_v1(monkeypatch, tmpdir):
    tmpdir.chdir()

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute1)
    charge, spin_multiplicity = check_charge_and_spin(TEST_ATOMS)
    output = relax_job(
        TEST_ATOMS,
        charge,
        spin_multiplicity,
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
    )

    assert output["atoms"] != TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826330655069403)

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(
        os.path.join(QCHEM_DIR, "mol.qin.basic.sella_opt_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_relax_job_v2(monkeypatch, tmpdir):
    tmpdir.chdir()
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute2)
    charge, spin_multiplicity = check_charge_and_spin(TEST_ATOMS, charge=-1)
    output = relax_job(
        TEST_ATOMS,
        charge,
        spin_multiplicity,
        method="b97mv",
        pcm_dielectric="3.0",
        opt_swaps={"max_steps": 1},
    )

    assert output["atoms"] != TEST_ATOMS
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
        os.path.join(QCHEM_DIR, "mol.qin.intermediate.sella_opt_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_relax_job_v3(monkeypatch, tmpdir):
    tmpdir.chdir()
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute3)
    overwrite_inputs = {"rem": {"mem_total": "170000"}}
    charge, spin_multiplicity = check_charge_and_spin(TEST_ATOMS)
    output = relax_job(
        TEST_ATOMS,
        charge,
        spin_multiplicity,
        scf_algorithm="gdm",
        overwrite_inputs=overwrite_inputs,
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
    )

    assert output["atoms"] != TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826311086011256)


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_relax_job_v4(tmpdir):
    tmpdir.chdir()
    with pytest.raises(ValueError):
        relax_job(TEST_ATOMS, 0, 1, pcm_dielectric="3.0", smd_solvent="water")


def test_freq_job_v1(monkeypatch, tmpdir):
    tmpdir.chdir()
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute5)
    charge, spin_multiplicity = check_charge_and_spin(TEST_ATOMS, charge=-1)
    output = freq_job(
        TEST_ATOMS,
        charge,
        spin_multiplicity,
        scf_algorithm="diis",
        method="b97mv",
        basis="def2-svpd",
    )

    assert output["atoms"] == TEST_ATOMS
    assert output["charge"] == -1
    assert output["spin_multiplicity"] == 2
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 77
    assert output["parameters"]["charge"] == -1
    assert output["parameters"]["spin_multiplicity"] == 2
    assert output["results"]["energy"] == pytest.approx(-605.6859554019 * units.Hartree)
    assert output["results"]["hessian"] is not None
    assert output["results"]["qc_output"]["enthalpy"] == pytest.approx(
        61.047 * (units.kcal / units.mol)
    )


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_ts_job_v1(monkeypatch, tmpdir):
    tmpdir.chdir()

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute1)
    charge, spin_multiplicity = check_charge_and_spin(TEST_ATOMS)
    output = ts_job(
        TEST_ATOMS,
        charge,
        spin_multiplicity,
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
    )

    assert output["atoms"] != TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826330655069403)

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(
        os.path.join(QCHEM_DIR, "mol.qin.basic.sella_TSopt_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_ts_job_v2(monkeypatch, tmpdir):
    tmpdir.chdir()
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute2)
    charge, spin_multiplicity = check_charge_and_spin(TEST_ATOMS, charge=-1)
    output = ts_job(
        TEST_ATOMS,
        charge,
        spin_multiplicity,
        method="b97mv",
        pcm_dielectric="3.0",
        opt_swaps={"max_steps": 1},
    )

    assert output["atoms"] != TEST_ATOMS
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
        os.path.join(QCHEM_DIR, "mol.qin.intermediate.sella_TSopt_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_ts_job_v3(monkeypatch, tmpdir):
    tmpdir.chdir()
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute3)
    overwrite_inputs = {"rem": {"mem_total": "170000"}}
    charge, spin_multiplicity = check_charge_and_spin(TEST_ATOMS)
    output = ts_job(
        TEST_ATOMS,
        charge,
        spin_multiplicity,
        scf_algorithm="gdm",
        overwrite_inputs=overwrite_inputs,
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
    )

    assert output["atoms"] != TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826311086011256)


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_ts_job_v4(tmpdir):
    tmpdir.chdir()
    with pytest.raises(ValueError):
        ts_job(TEST_ATOMS, 0, 1, pcm_dielectric="3.0", smd_solvent="water")

    with pytest.raises(ValueError):
        ts_job(
            TEST_ATOMS,
            0,
            1,
            pcm_dielectric="3.0",
            smd_solvent="water",
            opt_swaps={"optimizer": FIRE},
        )


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_irc_job_v1(monkeypatch, tmpdir):
    tmpdir.chdir()

    monkeypatch.setattr(QChem, "read_results", mock_read)
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute4)

    charge, spin_multiplicity = check_charge_and_spin(TEST_ATOMS)
    output = irc_job(
        TEST_ATOMS,
        charge,
        spin_multiplicity,
        "forward",
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
    )

    assert output["atoms"] != TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(
        os.path.join(QCHEM_DIR, "mol.qin.basic.sella_IRC_forward_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)

    charge, spin_multiplicity = check_charge_and_spin(TEST_ATOMS)
    output = irc_job(
        TEST_ATOMS,
        charge,
        spin_multiplicity,
        "reverse",
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
    )

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(
        os.path.join(QCHEM_DIR, "mol.qin.basic.sella_IRC_reverse_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)

    overwrite_inputs = {"rem": {"mem_total": "170000"}}
    output = irc_job(
        TEST_ATOMS,
        charge,
        spin_multiplicity,
        "reverse",
        scf_algorithm="gdm",
        overwrite_inputs=overwrite_inputs,
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
    )

    assert output["atoms"] != TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_irc_job_v2(tmpdir):
    tmpdir.chdir()
    with pytest.raises(ValueError):
        irc_job(TEST_ATOMS, 0, 1, "straight")

    with pytest.raises(ValueError):
        irc_job(
            TEST_ATOMS,
            0,
            1,
            "forward",
            pcm_dielectric="3.0",
            smd_solvent="water",
        )

    with pytest.raises(ValueError):
        irc_job(
            TEST_ATOMS,
            0,
            1,
            "forward",
            pcm_dielectric="3.0",
            smd_solvent="water",
            opt_swaps={"optimizer": FIRE},
        )


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_quasi_irc_job(monkeypatch, tmpdir):
    tmpdir.chdir()

    monkeypatch.setattr(QChem, "read_results", mock_read)
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute4)

    shared_kwargs = {"basis": "def2-tzvpd"}
    relax_opt_swaps = {"max_steps": 5}

    charge, spin_multiplicity = check_charge_and_spin(TEST_ATOMS)
    output = quasi_irc_job(
        TEST_ATOMS,
        charge,
        spin_multiplicity,
        "forward",
        shared_kwargs=shared_kwargs,
        relax_opt_swaps=relax_opt_swaps,
    )

    assert output["atoms"] != TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(
        os.path.join(QCHEM_DIR, "mol.qin.basic.quasi_irc_forward")
    )
    qcinput_nearly_equal(qcin, ref_qcin)

    shared_kwargs = {
        "basis": "def2-svpd",
        "scf_algorithm": "gdm",
    }
    irc_opt_swaps = {"max_steps": 6}
    relax_opt_swaps = {"max_steps": 6}

    output = quasi_irc_job(
        TEST_ATOMS,
        -1,
        2,
        "reverse",
        shared_kwargs=shared_kwargs,
        irc_opt_swaps=irc_opt_swaps,
        relax_opt_swaps=relax_opt_swaps,
    )

    assert output["atoms"] != TEST_ATOMS
    assert output["charge"] == -1
    assert output["spin_multiplicity"] == 2
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 77
    assert output["parameters"]["charge"] == -1
    assert output["parameters"]["spin_multiplicity"] == 2

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(os.path.join(QCHEM_DIR, "mol.qin.quasi_irc_reverse"))
    qcinput_nearly_equal(qcin, ref_qcin)


def test_internal_relax_job(monkeypatch, tmpdir):
    tmpdir.chdir()

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute1)
    charge, spin_multiplicity = check_charge_and_spin(TEST_ATOMS)
    output = internal_relax_job(TEST_ATOMS, charge, spin_multiplicity)
    assert output["atoms"] == TEST_ATOMS
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1
    assert output["results"]["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert output["results"]["forces"][0][0] == pytest.approx(-1.3826330655069403)

    qcin = QCInput.from_file("mol.qin.gz").as_dict()
    assert qcin["rem"]["basis"] == "def2-svpd"
    assert qcin["rem"]["geom_opt_max_cycles"] == "200"
    assert qcin["rem"]["job_type"] == "opt"
