import os
from pathlib import Path

import pytest
from ase import units
from ase.io import read
from monty.io import zopen
from pymatgen.io.qchem.inputs import QCInput

from quacc.calculators.qchem import QChem

FILE_DIR = Path(__file__).resolve().parent


@pytest.fixture()
def test_atoms():
    return read(FILE_DIR / "test.xyz")


@pytest.fixture()
def os_atoms():
    return read(FILE_DIR / "OS_test.xyz")


def test_qchem_write_input_basic(tmpdir, test_atoms):
    tmpdir.chdir()
    calc = QChem(
        test_atoms, rem={"method": "wb97xv", "basis": "def2-svp", "scf_guess": None}
    )
    assert calc.parameters["charge"] == 0
    assert calc.parameters["spin_multiplicity"] == 1
    assert calc.parameters["rem"] == {"method": "wb97xv", "basis": "def2-svp"}
    calc.write_input(test_atoms)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(str(FILE_DIR / "examples" / "basic" / "mol.qin"))
    assert qcinp.as_dict() == ref_qcinp.as_dict()
    assert not Path(FILE_DIR / "53.0").exists()

    with pytest.raises(NotImplementedError):
        QChem(
            test_atoms,
            rem={"method": "wb97xv", "basis": "def2-svp"},
            directory="notsupported",
        )


# def test_qchem_write_input_intermediate(tmpdir, test_atoms):
#     tmpdir.chdir()
#     params = {"dft_rung": 3, "pcm_dielectric": "3.0"}
#     calc = QChem(
#         test_atoms,
#         basis_set="def2-svpd",
#         cores=40,
#         charge=-1,
#         spin_multiplicity=2,
#         qchem_input_params=params,
#     )
#     assert calc.parameters["cores"] == 40
#     assert calc.parameters["charge"] == -1
#     assert calc.parameters["spin_multiplicity"] == 2
#     assert calc.parameters["dft_rung"] == 3
#     assert calc.parameters["basis_set"] == "def2-svpd"
#     assert calc.parameters["pcm_dielectric"] == "3.0"
#     assert calc.parameters["scf_algorithm"] == "diis"
#     calc.write_input(test_atoms)
#     qcinp = QCInput.from_file("mol.qin")
#     ref_qcinp = QCInput.from_file(
#         str(FILE_DIR / "examples" / "intermediate" / "mol.qin")
#     )
#     assert qcinp.as_dict() == ref_qcinp.as_dict()


# def test_qchem_write_input_advanced(tmpdir, test_atoms):
#     tmpdir.chdir()
#     params = {
#         "smd_solvent": "water",
#         "overwrite_inputs": {"rem": {"method": "b97mv", "mem_total": "170000"}},
#     }
#     calc = QChem(
#         test_atoms,
#         basis_set="def2-svpd",
#         scf_algorithm="gdm",
#         cores=40,
#         charge=-1,
#         spin_multiplicity=2,
#         qchem_input_params=params,
#     )
#     assert calc.parameters["cores"] == 40
#     assert calc.parameters["charge"] == -1
#     assert calc.parameters["spin_multiplicity"] == 2
#     assert calc.parameters["scf_algorithm"] == "gdm"
#     assert calc.parameters["basis_set"] == "def2-svpd"
#     assert calc.parameters["smd_solvent"] == "water"
#     assert calc.parameters["overwrite_rem_method"] == "b97mv"
#     assert calc.parameters["overwrite_rem_mem_total"] == "170000"
#     assert "method" not in calc.parameters
#     calc.write_input(test_atoms)
#     qcinp = QCInput.from_file("mol.qin")
#     ref_qcinp = QCInput.from_file(str(FILE_DIR / "examples" / "advanced" / "mol.qin"))
#     assert qcinp.as_dict() == ref_qcinp.as_dict()


# def test_qchem_write_input_open_shell_and_different_charges(tmpdir, os_atoms):
#     tmpdir.chdir()
#     calc = QChem(os_atoms, spin_multiplicity=2, cores=40)
#     assert calc.parameters["cores"] == 40
#     assert calc.parameters["charge"] == 0
#     assert calc.parameters["spin_multiplicity"] == 2
#     calc.write_input(os_atoms)
#     qcinp = QCInput.from_file("mol.qin")
#     ref_qcinp = QCInput.from_file(str(FILE_DIR / "examples" / "OSDC1.qin"))
#     assert qcinp.as_dict() == ref_qcinp.as_dict()

#     calc = QChem(os_atoms, cores=40, charge=0, spin_multiplicity=4)
#     assert calc.parameters["cores"] == 40
#     assert calc.parameters["charge"] == 0
#     assert calc.parameters["spin_multiplicity"] == 4
#     calc.write_input(os_atoms)
#     qcinp = QCInput.from_file("mol.qin")
#     ref_qcinp = QCInput.from_file(str(FILE_DIR / "examples" / "OSDC2.qin"))
#     assert qcinp.as_dict() == ref_qcinp.as_dict()

#     calc = QChem(os_atoms, cores=40, charge=1)
#     assert calc.parameters["cores"] == 40
#     assert calc.parameters["charge"] == 1
#     assert calc.parameters["spin_multiplicity"] == 1
#     calc.write_input(os_atoms)
#     qcinp = QCInput.from_file("mol.qin")
#     ref_qcinp = QCInput.from_file(str(FILE_DIR / "examples" / "OSDC3.qin"))
#     assert qcinp.as_dict() == ref_qcinp.as_dict()

#     calc = QChem(os_atoms, cores=40, charge=1, spin_multiplicity=3)
#     assert calc.parameters["cores"] == 40
#     assert calc.parameters["charge"] == 1
#     assert calc.parameters["spin_multiplicity"] == 3
#     calc.write_input(os_atoms)
#     qcinp = QCInput.from_file("mol.qin")
#     ref_qcinp = QCInput.from_file(str(FILE_DIR / "examples" / "OSDC4.qin"))
#     assert qcinp.as_dict() == ref_qcinp.as_dict()


# def test_qchem_write_input_freq(tmpdir, test_atoms):
#     tmpdir.chdir()
#     params = {"dft_rung": 3, "pcm_dielectric": "3.0"}
#     calc = QChem(
#         test_atoms,
#         job_type="freq",
#         basis_set="def2-svpd",
#         cores=40,
#         charge=-1,
#         spin_multiplicity=2,
#         qchem_input_params=params,
#     )
#     assert calc.parameters["cores"] == 40
#     assert calc.parameters["charge"] == -1
#     assert calc.parameters["spin_multiplicity"] == 2
#     assert calc.parameters["dft_rung"] == 3
#     assert calc.parameters["basis_set"] == "def2-svpd"
#     assert calc.parameters["pcm_dielectric"] == "3.0"
#     assert calc.parameters["scf_algorithm"] == "diis"
#     calc.write_input(test_atoms)
#     qcinp = QCInput.from_file("mol.qin")
#     ref_qcinp = QCInput.from_file(str(FILE_DIR / "examples" / "freq" / "mol.qin"))
#     assert qcinp.as_dict() == ref_qcinp.as_dict()


# def test_qchem_read_results_basic_and_write_53(tmpdir, test_atoms):
#     calc = QChem(test_atoms, cores=40)
#     os.chdir(FILE_DIR / "examples" / "basic")
#     calc.read_results()
#     tmpdir.chdir()

#     assert calc.results["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
#     assert calc.results["forces"][0][0] == pytest.approx(-1.3826330655069403)
#     assert calc._prev_orbital_coeffs is not None

#     calc.write_input(test_atoms)
#     assert Path(tmpdir, "53.0").exists()
#     with zopen("53.0", mode="rb") as new_file:
#         new_binary = new_file.read()
#         with zopen(
#             os.path.join(FILE_DIR, "examples", "basic", "53.0"), mode="rb"
#         ) as old_file:
#             old_binary = old_file.read()
#             assert new_binary == old_binary
#     qcinp = QCInput.from_file("mol.qin")
#     assert qcinp.rem.get("scf_guess") == "read"


# def test_qchem_read_results_intermediate(tmpdir, test_atoms):
#     tmpdir.chdir()
#     calc = QChem(test_atoms, cores=40)
#     os.chdir(FILE_DIR / "examples" / "intermediate")
#     calc.read_results()
#     tmpdir.chdir()

#     assert calc.results["energy"] == pytest.approx(-605.6859554025 * units.Hartree)
#     assert calc.results["forces"][0][0] == pytest.approx(-0.6955571014353796)
#     assert calc._prev_orbital_coeffs is not None


# def test_qchem_read_results_advanced(tmpdir, test_atoms):
#     tmpdir.chdir()
#     calc = QChem(test_atoms, cores=40)
#     os.chdir(FILE_DIR / "examples" / "advanced")
#     calc.read_results()
#     tmpdir.chdir()

#     assert calc.results["energy"] == pytest.approx(-605.7310332390 * units.Hartree)
#     assert calc.results["forces"][0][0] == pytest.approx(-0.4270884974249971)
#     assert calc._prev_orbital_coeffs is not None
#     assert calc.results.get("hessian") is None


# def test_qchem_read_results_freq(tmpdir, test_atoms):
#     calc = QChem(test_atoms, job_type="freq", cores=40)
#     os.chdir(FILE_DIR / "examples" / "freq")
#     calc.read_results()
#     tmpdir.chdir()

#     assert calc.results["energy"] == pytest.approx(-605.6859554025 * units.Hartree)
#     assert calc.results.get("forces") is None
#     assert calc._prev_orbital_coeffs is not None
#     assert len(calc.results["hessian"]) == 42
#     assert len(calc.results["hessian"][0]) == 42
#     assert calc.results["qc_output"]["frequencies"][0] == -340.2
#     assert len(calc.results["qc_output"]["frequencies"]) == 36
#     assert len(calc.results["qc_output"]["frequency_mode_vectors"]) == 36
#     assert len(calc.results["qc_output"]["frequency_mode_vectors"][0]) == 14
#     assert len(calc.results["qc_output"]["frequency_mode_vectors"][0][0]) == 3
