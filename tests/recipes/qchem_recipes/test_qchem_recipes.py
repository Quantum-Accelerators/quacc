import pytest

pytest.importorskip("sella")


@pytest.fixture()
def qchem_dir():
    from pathlib import Path

    FILE_DIR = Path(__file__).resolve().parent
    return FILE_DIR / "qchem_examples"


@pytest.fixture()
def test_atoms():
    from pathlib import Path

    from ase.io import read

    FILE_DIR = Path(__file__).resolve().parent

    return read(FILE_DIR / "test.xyz")


@pytest.fixture()
def os_atoms():
    from pathlib import Path

    from ase.io import read

    FILE_DIR = Path(__file__).resolve().parent

    return read(FILE_DIR / "OS_test.xyz")


def setup_module():
    from quacc import SETTINGS

    SETTINGS.CHECK_CONVERGENCE = False


def teardown_module():
    from quacc import SETTINGS

    SETTINGS.CHECK_CONVERGENCE = True


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
    from pathlib import Path
    from shutil import copy

    FILE_DIR = Path(__file__).resolve().parent
    qchem_dir = FILE_DIR / "qchem_examples"
    copy(qchem_dir / "mol.qout.basic", "mol.qout")
    copy(qchem_dir / "131.0.basic", "131.0")
    copy(qchem_dir / "53.0.basic", "53.0")
    copy(qchem_dir / "custodian.json", "custodian.json")


def mock_execute2(_self, **kwargs):
    from pathlib import Path
    from shutil import copy

    FILE_DIR = Path(__file__).resolve().parent
    qchem_dir = FILE_DIR / "qchem_examples"
    copy(qchem_dir / "mol.qout.intermediate", "mol.qout")
    copy(qchem_dir / "131.0.intermediate", "131.0")
    copy(qchem_dir / "53.0.intermediate", "53.0")
    copy(qchem_dir / "custodian.json", "custodian.json")


def mock_execute3(_self, **kwargs):
    from pathlib import Path
    from shutil import copy

    FILE_DIR = Path(__file__).resolve().parent
    qchem_dir = FILE_DIR / "qchem_examples"
    copy(qchem_dir / "mol.qout.alternate", "mol.qout")
    copy(qchem_dir / "131.0.alternate", "131.0")
    copy(qchem_dir / "53.0.alternate", "53.0")
    copy(qchem_dir / "custodian.json", "custodian.json")


def mock_execute4(self, **kwargs):
    from ase.calculators.lj import LennardJones
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.io.qchem.inputs import QCInput

    qcin = QCInput.from_file("mol.qin")
    mol = qcin.molecule
    atoms = AseAtomsAdaptor.get_atoms(mol)
    atoms.calc = LennardJones()
    atoms.get_potential_energy()
    self.results = atoms.calc.results


def mock_execute5(_self, **kwargs):
    from pathlib import Path
    from shutil import copy

    FILE_DIR = Path(__file__).resolve().parent
    qchem_dir = FILE_DIR / "qchem_examples"

    copy(qchem_dir / "mol.qout.freq", "mol.qout")
    copy(qchem_dir / "132.0.freq", "132.0")
    copy(qchem_dir / "53.0.freq", "53.0")
    copy(qchem_dir / "custodian.json", "custodian.json")


def mock_read(self, **kwargs):
    if self.results is None:
        raise RuntimeError("Results should not be None here.")


def test_static_job_v1(monkeypatch, tmpdir, qchem_dir, test_atoms):
    from ase import units
    from ase.calculators.calculator import FileIOCalculator
    from pymatgen.io.qchem.inputs import QCInput

    from quacc.atoms.core import check_charge_and_spin
    from quacc.recipes.qchem.core import static_job

    tmpdir.chdir()
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute1)
    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = static_job(test_atoms, charge, spin_multiplicity)
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
    ref_qcin = QCInput.from_file(str(qchem_dir / "mol.qin.basic"))
    qcinput_nearly_equal(qcin, ref_qcin)
    qcinput_nearly_equal(ref_qcin, QCInput.from_dict(output["results"]["qc_input"]))


def test_static_job_v2(monkeypatch, tmpdir, qchem_dir, test_atoms):
    from ase import units
    from ase.calculators.calculator import FileIOCalculator
    from pymatgen.io.qchem.inputs import QCInput

    from quacc.atoms.core import check_charge_and_spin
    from quacc.recipes.qchem.core import static_job

    tmpdir.chdir()

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute2)
    charge, spin_multiplicity = check_charge_and_spin(test_atoms, charge=-1)
    output = static_job(
        test_atoms,
        charge,
        spin_multiplicity,
        method="b97mv",
        basis="def2-svpd",
        pcm_dielectric="3.0",
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
    ref_qcin = QCInput.from_file(str(qchem_dir / "mol.qin.intermediate"))
    qcinput_nearly_equal(qcin, ref_qcin)
    qcinput_nearly_equal(ref_qcin, QCInput.from_dict(output["results"]["qc_input"]))


def test_static_job_v3(monkeypatch, tmpdir, qchem_dir, test_atoms):
    from ase import units
    from ase.calculators.calculator import FileIOCalculator
    from pymatgen.io.qchem.inputs import QCInput

    from quacc.atoms.core import check_charge_and_spin
    from quacc.recipes.qchem.core import static_job

    tmpdir.chdir()

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute3)
    overwrite_inputs = {"rem": {"mem_total": "170000"}}
    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = static_job(
        test_atoms,
        charge,
        spin_multiplicity,
        scf_algorithm="gdm",
        overwrite_inputs=overwrite_inputs,
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
    ref_qcin = QCInput.from_file(str(qchem_dir / "mol.qin.alternate"))
    qcinput_nearly_equal(qcin, ref_qcin)
    qcinput_nearly_equal(ref_qcin, QCInput.from_dict(output["results"]["qc_input"]))


def test_static_job_v4(monkeypatch, tmpdir, os_atoms):
    from ase.calculators.calculator import FileIOCalculator

    from quacc.atoms.core import check_charge_and_spin
    from quacc.calculators.qchem import QChem
    from quacc.recipes.qchem.core import static_job

    tmpdir.chdir()
    monkeypatch.setattr(QChem, "read_results", mock_read)
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute4)
    charge, spin_multiplicity = check_charge_and_spin(os_atoms)
    assert static_job(os_atoms, charge, spin_multiplicity)


def test_static_job_v5(tmpdir, test_atoms):
    from quacc.recipes.qchem.core import static_job

    tmpdir.chdir()

    with pytest.raises(ValueError):
        static_job(test_atoms, 0, 1, pcm_dielectric="3.0", smd_solvent="water")


def test_relax_job_v1(monkeypatch, tmpdir, qchem_dir, test_atoms):
    from ase import units
    from ase.calculators.calculator import FileIOCalculator
    from pymatgen.io.qchem.inputs import QCInput

    from quacc.atoms.core import check_charge_and_spin
    from quacc.recipes.qchem.core import relax_job

    tmpdir.chdir()

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute1)
    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = relax_job(
        test_atoms,
        charge,
        spin_multiplicity,
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
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
    ref_qcin = QCInput.from_file(str(qchem_dir / "mol.qin.basic.sella_opt_iter1"))
    qcinput_nearly_equal(qcin, ref_qcin)
    assert len(output["results"]["qc_input"]) > 1


def test_relax_job_v2(monkeypatch, tmpdir, qchem_dir, test_atoms):
    from ase import units
    from ase.calculators.calculator import FileIOCalculator
    from pymatgen.io.qchem.inputs import QCInput

    from quacc.atoms.core import check_charge_and_spin
    from quacc.recipes.qchem.core import relax_job

    tmpdir.chdir()
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute2)
    charge, spin_multiplicity = check_charge_and_spin(test_atoms, charge=-1)
    output = relax_job(
        test_atoms,
        charge,
        spin_multiplicity,
        method="b97mv",
        pcm_dielectric="3.0",
        opt_swaps={"max_steps": 1},
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
        str(qchem_dir / "mol.qin.intermediate.sella_opt_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)


def test_relax_job_v3(monkeypatch, tmpdir, test_atoms):
    from ase import units
    from ase.calculators.calculator import FileIOCalculator

    from quacc.atoms.core import check_charge_and_spin
    from quacc.recipes.qchem.core import relax_job

    tmpdir.chdir()
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute3)
    overwrite_inputs = {"rem": {"mem_total": "170000"}}
    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = relax_job(
        test_atoms,
        charge,
        spin_multiplicity,
        scf_algorithm="gdm",
        overwrite_inputs=overwrite_inputs,
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
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


def test_relax_job_v4(tmpdir, test_atoms):
    from quacc.recipes.qchem.core import relax_job

    tmpdir.chdir()
    with pytest.raises(ValueError):
        relax_job(test_atoms, 0, 1, pcm_dielectric="3.0", smd_solvent="water")


def test_freq_job_v1(monkeypatch, tmpdir, test_atoms):
    from ase import units
    from ase.calculators.calculator import FileIOCalculator

    from quacc.atoms.core import check_charge_and_spin
    from quacc.recipes.qchem.core import freq_job

    tmpdir.chdir()
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute5)
    charge, spin_multiplicity = check_charge_and_spin(test_atoms, charge=-1)
    output = freq_job(
        test_atoms,
        charge,
        spin_multiplicity,
        scf_algorithm="diis",
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
    assert output["results"]["hessian"] is not None
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


def test_ts_job_v1(monkeypatch, tmpdir, qchem_dir, test_atoms):
    from ase import units
    from ase.calculators.calculator import FileIOCalculator
    from pymatgen.io.qchem.inputs import QCInput

    from quacc.atoms.core import check_charge_and_spin
    from quacc.recipes.qchem.ts import ts_job

    tmpdir.chdir()

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute1)
    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = ts_job(
        test_atoms,
        charge,
        spin_multiplicity,
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
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
    ref_qcin = QCInput.from_file(str(qchem_dir / "mol.qin.basic.sella_TSopt_iter1"))
    qcinput_nearly_equal(qcin, ref_qcin)


def test_ts_job_v2(monkeypatch, tmpdir, qchem_dir, test_atoms):
    from ase import units
    from ase.calculators.calculator import FileIOCalculator
    from pymatgen.io.qchem.inputs import QCInput

    from quacc.atoms.core import check_charge_and_spin
    from quacc.recipes.qchem.ts import ts_job

    tmpdir.chdir()
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute2)
    charge, spin_multiplicity = check_charge_and_spin(test_atoms, charge=-1)
    output = ts_job(
        test_atoms,
        charge,
        spin_multiplicity,
        method="b97mv",
        pcm_dielectric="3.0",
        opt_swaps={"max_steps": 1},
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
        str(qchem_dir / "mol.qin.intermediate.sella_TSopt_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)


def test_ts_job_v3(monkeypatch, tmpdir, test_atoms):
    from ase import units
    from ase.calculators.calculator import FileIOCalculator

    from quacc.atoms.core import check_charge_and_spin
    from quacc.recipes.qchem.ts import ts_job

    tmpdir.chdir()
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute3)
    overwrite_inputs = {"rem": {"mem_total": "170000"}}
    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = ts_job(
        test_atoms,
        charge,
        spin_multiplicity,
        scf_algorithm="gdm",
        overwrite_inputs=overwrite_inputs,
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
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


def test_ts_job_v4(tmpdir, test_atoms):
    from ase.optimize import FIRE

    from quacc.recipes.qchem.ts import ts_job

    tmpdir.chdir()
    with pytest.raises(ValueError):
        ts_job(test_atoms, 0, 1, pcm_dielectric="3.0", smd_solvent="water")

    with pytest.raises(ValueError):
        ts_job(
            test_atoms,
            0,
            1,
            pcm_dielectric="3.0",
            smd_solvent="water",
            opt_swaps={"optimizer": FIRE},
        )


def test_irc_job_v1(monkeypatch, tmpdir, qchem_dir, test_atoms):
    from ase.calculators.calculator import FileIOCalculator
    from pymatgen.io.qchem.inputs import QCInput

    from quacc.atoms.core import check_charge_and_spin
    from quacc.calculators.qchem import QChem
    from quacc.recipes.qchem.ts import irc_job

    tmpdir.chdir()

    monkeypatch.setattr(QChem, "read_results", mock_read)
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute4)

    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = irc_job(
        test_atoms,
        charge,
        spin_multiplicity,
        "forward",
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
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
        str(qchem_dir / "mol.qin.basic.sella_IRC_forward_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)

    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = irc_job(
        test_atoms,
        charge,
        spin_multiplicity,
        "reverse",
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
    )

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(
        str(qchem_dir / "mol.qin.basic.sella_IRC_reverse_iter1")
    )
    qcinput_nearly_equal(qcin, ref_qcin)

    overwrite_inputs = {"rem": {"mem_total": "170000"}}
    output = irc_job(
        test_atoms,
        charge,
        spin_multiplicity,
        "reverse",
        scf_algorithm="gdm",
        overwrite_inputs=overwrite_inputs,
        basis="def2-tzvpd",
        opt_swaps={"max_steps": 1},
    )

    assert output["atoms"] != test_atoms
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1


def test_irc_job_v2(tmpdir, test_atoms):
    from ase.optimize import FIRE

    from quacc.recipes.qchem.ts import irc_job

    tmpdir.chdir()
    with pytest.raises(ValueError):
        irc_job(test_atoms, 0, 1, "straight")

    with pytest.raises(ValueError):
        irc_job(
            test_atoms,
            0,
            1,
            "forward",
            pcm_dielectric="3.0",
            smd_solvent="water",
        )

    with pytest.raises(ValueError):
        irc_job(
            test_atoms,
            0,
            1,
            "forward",
            pcm_dielectric="3.0",
            smd_solvent="water",
            opt_swaps={"optimizer": FIRE},
        )


def test_quasi_irc_job(monkeypatch, tmpdir, qchem_dir, test_atoms):
    from ase.calculators.calculator import FileIOCalculator
    from pymatgen.io.qchem.inputs import QCInput

    from quacc.atoms.core import check_charge_and_spin
    from quacc.calculators.qchem import QChem
    from quacc.recipes.qchem.ts import quasi_irc_job

    tmpdir.chdir()

    monkeypatch.setattr(QChem, "read_results", mock_read)
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute4)

    relax_opt_swaps = {"max_steps": 5}

    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = quasi_irc_job(
        test_atoms,
        charge,
        spin_multiplicity,
        "forward",
        basis="def2-tzvpd",
        relax_opt_swaps=relax_opt_swaps,
    )

    assert output["atoms"] != test_atoms
    assert output["charge"] == 0
    assert output["spin_multiplicity"] == 1
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 76
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 1

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(str(qchem_dir / "mol.qin.basic.quasi_irc_forward"))
    qcinput_nearly_equal(qcin, ref_qcin)

    irc_opt_swaps = {"max_steps": 6}
    relax_opt_swaps = {"max_steps": 6}

    output = quasi_irc_job(
        test_atoms,
        -1,
        2,
        "reverse",
        basis="def2-svpd",
        scf_algorithm="gdm",
        irc_opt_swaps=irc_opt_swaps,
        relax_opt_swaps=relax_opt_swaps,
    )

    assert output["atoms"] != test_atoms
    assert output["charge"] == -1
    assert output["spin_multiplicity"] == 2
    assert output["formula_alphabetical"] == "C4 H4 O6"
    assert output["nelectrons"] == 77
    assert output["parameters"]["charge"] == -1
    assert output["parameters"]["spin_multiplicity"] == 2

    qcin = QCInput.from_file("mol.qin.gz")
    ref_qcin = QCInput.from_file(str(qchem_dir / "mol.qin.quasi_irc_reverse"))
    qcinput_nearly_equal(qcin, ref_qcin)


def test_internal_relax_job(monkeypatch, tmpdir, test_atoms):
    from ase import units
    from ase.calculators.calculator import FileIOCalculator
    from pymatgen.io.qchem.inputs import QCInput

    from quacc.atoms.core import check_charge_and_spin
    from quacc.recipes.qchem.core import internal_relax_job

    tmpdir.chdir()

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute1)
    charge, spin_multiplicity = check_charge_and_spin(test_atoms)
    output = internal_relax_job(test_atoms, charge, spin_multiplicity)
    assert output["atoms"] == test_atoms
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
