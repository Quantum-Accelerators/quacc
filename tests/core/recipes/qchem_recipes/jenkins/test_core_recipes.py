from pathlib import Path

import pytest
from ase import units
from ase.io import read
from pymatgen.io.qchem.inputs import QCInput

from quacc import SETTINGS
from quacc.atoms.core import check_charge_and_spin
from quacc.recipes.qchem.core import static_job

FILE_DIR = Path(__file__).parent
QCHEM_DIR = FILE_DIR / "qchem_examples"

DEFAULT_SETTINGS = SETTINGS.model_copy()


@pytest.fixture()
def test_atoms():
    return read(FILE_DIR / "test.xyz")


@pytest.fixture()
def os_atoms():
    return read(FILE_DIR / "OS_test.xyz")


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


def test_static_job_v1(monkeypatch, tmp_path, test_atoms):
    monkeypatch.chdir(tmp_path)
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
