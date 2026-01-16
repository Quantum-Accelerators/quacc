from __future__ import annotations

from logging import WARNING, getLogger
from shutil import which

import pytest

from quacc import get_settings

settings = get_settings()
pytestmark = pytest.mark.skipif(
    which(str(settings.ESPRESSO_BIN_DIR / settings.ESPRESSO_BINARIES["pw"])) is None
    or which(str(settings.ESPRESSO_BIN_DIR / settings.ESPRESSO_BINARIES["ph"])) is None,
    reason="QE not installed",
)

from pathlib import Path

from ase.build import bulk
from ase.io.espresso import read_fortran_namelist
from monty.io import zopen
from monty.shutil import decompress_file
from numpy.testing import assert_allclose, assert_array_equal

from quacc import change_settings
from quacc.recipes.espresso.core import non_scf_job, relax_job, static_job
from quacc.recipes.espresso.phonons import (
    dvscf_q2r_job,
    matdyn_job,
    phonon_dos_flow,
    phonon_job,
    postahc_job,
    q2r_job,
)
from quacc.utils.files import copy_decompress_files
from quacc.wflow_tools.job_argument import Copy

DATA_DIR = Path(__file__).parent / "data"

LOGGER = getLogger(__name__)
LOGGER.propagate = True


def test_phonon_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    with change_settings({"ESPRESSO_PSEUDO": tmp_path}):
        atoms = bulk("Li")

        copy_decompress_files(DATA_DIR, ["Li.upf.gz"], tmp_path)

        input_data = {
            "control": {"calculation": "scf"},
            "system": {"occupations": "smearing", "smearing": "cold", "degauss": 0.02},
            "electrons": {"mixing_mode": "TF", "mixing_beta": 0.7, "conv_thr": 1.0e-6},
        }

        ph_loose = {"inputph": {"tr2_ph": 1e-8}}

        pseudopotentials = {"Li": "Li.upf"}

        pw_results = static_job(
            atoms,
            input_data=input_data,
            pseudopotentials=pseudopotentials,
            kspacing=0.5,
        )

        ph_results = phonon_job(pw_results["dir_name"], input_data=ph_loose)

        assert_allclose(
            ph_results["results"][1]["atoms"].get_positions(),
            atoms.get_positions(),
            atol=1.0e-4,
        )
        # ph.x cell param are not defined to a very high level of accuracy,
        # atol = 1.0e-3 is needed here...
        assert_allclose(
            ph_results["results"][1]["atoms"].get_cell(), atoms.get_cell(), atol=1.0e-3
        )
        assert_array_equal(
            ph_results["results"][1]["atoms"].get_chemical_symbols(),
            atoms.get_chemical_symbols(),
        )

        sections = [
            "atoms",
            "eqpoints",
            "freqs",
            "kpoints",
            "mode_symmetries",
            "representations",
        ]

        for key in sections:
            assert key in ph_results["results"][1]

        ph_loose["inputph"]["recover"] = True

        recover_ph_results = phonon_job(ph_results["dir_name"], input_data=ph_loose)

        with zopen(Path(recover_ph_results["dir_name"], "ph.out.gz")) as f:
            lines = str(f.read())

        assert "Reading collected, re-writing distributed wavefunctions in" in lines
        assert "Restart after Phonon calculation" in lines
        assert "Representation     1      1 modes -  Done" in lines
        assert "Representation     2      1 modes -  Done" in lines
        assert "Representation     3      1 modes -  Done" in lines


def test_phonon_job_lqdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    with change_settings({"ESPRESSO_PSEUDO": tmp_path}):
        atoms = bulk("Li")

        copy_decompress_files(DATA_DIR, ["Li.upf.gz"], tmp_path)

        input_data = {
            "control": {"calculation": "scf"},
            "system": {"occupations": "smearing", "smearing": "cold", "degauss": 0.02},
            "electrons": {"mixing_mode": "TF", "mixing_beta": 0.7, "conv_thr": 1.0e-6},
        }

        ph_loose = {
            "inputph": {"tr2_ph": 1e-8, "lqdir": True, "qplot": True, "ldisp": True}
        }

        pseudopotentials = {"Li": "Li.upf"}

        pw_results = static_job(
            atoms,
            input_data=input_data,
            pseudopotentials=pseudopotentials,
            kspacing=0.5,
        )

        ph_results = phonon_job(
            pw_results["dir_name"],
            input_data=ph_loose,
            qpts=[(0, 0, 0, 1), (0.1, 0, 0, 1)],
        )

        assert_allclose(
            ph_results["results"][1]["atoms"].get_positions(),
            atoms.get_positions(),
            atol=1.0e-4,
        )
        # ph.x cell param are not defined to a very high level of accuracy,
        # atol = 1.0e-3 is needed here...
        assert_allclose(
            ph_results["results"][1]["atoms"].get_cell(), atoms.get_cell(), atol=1.0e-3
        )
        assert_array_equal(
            ph_results["results"][1]["atoms"].get_chemical_symbols(),
            atoms.get_chemical_symbols(),
        )

        sections = [
            "atoms",
            "eqpoints",
            "freqs",
            "kpoints",
            "mode_symmetries",
            "representations",
        ]

        for key in sections:
            assert key in ph_results["results"][1]

        ph_loose["inputph"]["recover"] = True

        recover_ph_results = phonon_job(
            ph_results["dir_name"],
            input_data=ph_loose,
            qpts=[(0, 0, 0, 1), (0.1, 0, 0, 1)],
        )

        with zopen(Path(recover_ph_results["dir_name"], "ph.out.gz")) as f:
            lines = str(f.read())

        assert "Reading collected, re-writing distributed wavefunctions in" in lines
        assert "Restart after Phonon calculation" in lines
        assert "Representation     1      3 modes -  Done" in lines
        assert "Representation     1      1 modes -  Done" in lines
        assert "Representation     2      2 modes -  Done" in lines


def test_phonon_job_list_to_do(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    with change_settings({"ESPRESSO_PSEUDO": tmp_path}):
        atoms = bulk("Li")

        copy_decompress_files(DATA_DIR, ["Li.upf.gz"], tmp_path)

        input_data = {
            "control": {"calculation": "scf"},
            "system": {"occupations": "smearing", "smearing": "cold", "degauss": 0.02},
            "electrons": {"mixing_mode": "TF", "mixing_beta": 0.7, "conv_thr": 1.0e-6},
        }
        ph_loose = {
            "inputph": {"tr2_ph": 1e-8, "qplot": True, "nat_todo": 1, "ldisp": True}
        }
        pseudopotentials = {"Li": "Li.upf"}

        pw_results = static_job(
            atoms,
            input_data=input_data,
            pseudopotentials=pseudopotentials,
            kspacing=0.5,
        )

        qpts = [(0, 0, 0, 1), (1 / 3, 0, 0, 1), (1 / 2, 0, 0, 1)]

        nat_todo = [1]

        ph_results = phonon_job(
            pw_results["dir_name"],
            input_data=ph_loose,
            qpts=qpts,
            nat_todo_indices=nat_todo,
        )

        assert_allclose(
            ph_results["results"][1]["atoms"].get_positions(),
            atoms.get_positions(),
            atol=1.0e-4,
        )
        # ph.x cell param are not defined to a very high level of accuracy,
        # atol = 1.0e-3 is needed here...
        assert_allclose(
            ph_results["results"][1]["atoms"].get_cell(), atoms.get_cell(), atol=1.0e-3
        )
        assert_array_equal(
            ph_results["results"][1]["atoms"].get_chemical_symbols(),
            atoms.get_chemical_symbols(),
        )

        sections = [
            "atoms",
            "eqpoints",
            "freqs",
            "kpoints",
            "mode_symmetries",
            "representations",
        ]

        for key in sections:
            assert key in ph_results["results"][1]


def test_q2r_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    with change_settings({"ESPRESSO_PSEUDO": tmp_path}):
        copy_decompress_files(DATA_DIR / "q2r_test", "matdyn", tmp_path)

        additional_cards = ["1 1 1", "1", "matdyn"]

        q2r_results = q2r_job(tmp_path, additional_cards=additional_cards)

        assert Path(q2r_results["dir_name"], "q2r.fc.gz").exists()

        decompress_file(Path(q2r_results["dir_name"], "q2r.in.gz"))

        with Path(q2r_results["dir_name"], "q2r.in").open() as f:
            recycled_input = read_fortran_namelist(f)

        assert recycled_input[0]["input"].pop("flfrc") == "q2r.fc"

        assert recycled_input[0]["input"] == {"fildyn": "matdyn"}
        assert recycled_input[1] == [*additional_cards, "EOF"]


def test_matdyn_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    with change_settings({"ESPRESSO_PSEUDO": tmp_path}):
        copy_decompress_files(DATA_DIR / "matdyn_test", "q2r.fc", tmp_path)

        input_data = {"input": {"dos": True, "nk1": 4, "nk2": 4, "nk3": 4}}
        matdyn_results = matdyn_job(tmp_path, input_data=input_data)

        assert Path(matdyn_results["dir_name"], "q2r.fc.gz").exists()
        assert Path(matdyn_results["dir_name"], "matdyn.dos.gz").exists()
        assert Path(matdyn_results["dir_name"], "matdyn.freq.gz").exists()
        assert Path(matdyn_results["dir_name"], "matdyn.modes.gz").exists()
        assert matdyn_results["results"]["matdyn_results"]["phonon_dos"].shape == (
            561,
            3,
        )

        decompress_file(Path(matdyn_results["dir_name"], "matdyn.in.gz"))

        with Path(matdyn_results["dir_name"], "matdyn.in").open() as f:
            recycled_input = read_fortran_namelist(f)

        assert recycled_input[0]["input"].pop("flfrc") == "q2r.fc"

        assert recycled_input[0]["input"]["fldos"] == "matdyn.dos"
        assert recycled_input[0]["input"]["flfrq"] == "matdyn.freq"
        assert recycled_input[0]["input"]["flvec"] == "matdyn.modes"
        assert recycled_input[0]["input"]["fleig"] == "matdyn.eig"


def test_phonon_dos_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    with change_settings({"ESPRESSO_PSEUDO": tmp_path}):
        copy_decompress_files(DATA_DIR, ["Li.upf.gz"], tmp_path)

        atoms = bulk("Li")
        input_data = {
            "control": {"calculation": "scf"},
            "system": {"occupations": "smearing", "smearing": "cold", "degauss": 0.02},
            "electrons": {"mixing_mode": "TF", "mixing_beta": 0.7, "conv_thr": 1.0e-6},
        }

        pseudopotentials = {"Li": "Li.upf"}

        relax_output = relax_job(
            atoms,
            input_data=input_data,
            pseudopotentials=pseudopotentials,
            kspacing=1.0,
        )
        assert phonon_dos_flow(prev_outdir=relax_output["dir_name"])


def test_phonon_calculation_spin_orbit_example_06(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    with change_settings({"ESPRESSO_PSEUDO": tmp_path}):
        # Create input structure
        pt_atoms = bulk("Pt")

        copy_decompress_files(DATA_DIR, ["Pt.rel-pz-n-rrkjus.UPF.gz"], tmp_path)

        # Set up pseudopotential
        pseudopotential = {"Pt": "Pt.rel-pz-n-rrkjus.UPF"}

        # Pt calculations
        pt_relax_params = {
            "input_data": {
                "cONTROL": {
                    "caLCULATION": "scf",
                    "RESTART_mode": "from_scratch",
                    "tPRNFOR": True,
                    "TstRESS": True,
                },
                "sYSTEM": {
                    "LSPINOrb": True,
                    "noNCOLIN": True,
                    "STARTING_MAGNETIZATION": 0.0,
                    "occupaTIONS": "smearing",
                    "DEGAUss": 0.02,
                    "SMEaring": "mv",
                    "eCUTWFC": 30.0,
                    "EcuTRHO": 250.0,
                },
                "ELECTrons": {"mixing_beta": 0.7, "conv_thr": 1.0e-8},
            },
            "pseudopotentials": pseudopotential,
            "kpts": (2, 2, 2),
        }
        pt_relax_results = relax_job(pt_atoms, **pt_relax_params)

        pt_phonon_x_params = {
            "input_data": {
                "inputph": {"aMass(1)": 195.078, "fildYn": "ptdyn", "TR2_ph": 1.0e-10}
            }
        }
        pt_phonon_x_results = phonon_job(
            pt_relax_results["dir_name"], **pt_phonon_x_params, qpts=(1.0, 0.0, 0.0)
        )

        with zopen(Path(pt_phonon_x_results["dir_name"], "ph.out.gz")) as f:
            lines = str(f.read())

        assert "Reading collected, re-writing distributed wavefunctions in" in lines
        assert "Non magnetic calculation with spin-orbit" in lines


def test_phonon_calculation_si_spin_orbit(tmp_path, monkeypatch, caplog):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    with change_settings({"ESPRESSO_PSEUDO": tmp_path}):
        copy_decompress_files(DATA_DIR, ["Si.rel-pbe.rrkj.UPF.gz"], tmp_path)

        si_atoms = bulk("Si")

        pseudopotential = {"Si": "Si.rel-pbe.rrkj.UPF"}

        si_relax_params = {
            "input_data": {
                "control": {"calculation": "scf", "restart_mode": "from_scratch"},
                "system": {
                    "ecutwfc": 20.0,
                    "noncolin": True,
                    "lspinorb": True,
                    "occupations": "fixed",
                },
                "electrons": {"mixing_beta": 0.7, "conv_thr": 1.0e-8},
            },
            "pseudopotentials": pseudopotential,
            "kpts": (2, 2, 2),
        }

        with caplog.at_level(WARNING):
            si_relax_results = relax_job(si_atoms, **si_relax_params)

            assert "The occupations are set to 'fixed'" in caplog.text
            caplog.clear()

        si_phonon_params = {
            "input_data": {
                "inputph": {
                    "tr2_ph": 1.0e-10,
                    "epsil": True,
                    "fildyn": "Sig.dyn",
                    "amass(1)": 28.0855,
                }
            }
        }
        si_phonon_results = phonon_job(
            si_relax_results["dir_name"], **si_phonon_params, qpts=(0.0, 0.0, 0.0)
        )

        assert (
            si_phonon_results["parameters"]["input_data"]["inputph"]["prefix"]
            == "pwscf"
        )

        assert (
            si_phonon_results["parameters"]["input_data"]["inputph"]["fildyn"]
            == "matdyn"
        )

        with zopen(Path(si_phonon_results["dir_name"], "ph.out.gz")) as f:
            lines = str(f.read())

        assert Path(
            si_phonon_results["dir_name"], "_ph0", "pwscf.phsave", "tensors.xml.gz"
        ).exists()

        assert "Reading collected, re-writing distributed wavefunctions in" in lines
        assert "Non magnetic calculation with spin-orbit" in lines


def test_phonon_induced_renormalization(tmp_path, monkeypatch, caplog):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    with change_settings({"ESPRESSO_PSEUDO": tmp_path}):
        copy_decompress_files(DATA_DIR, ["C.UPF.gz"], tmp_path)

        atoms = bulk("C")

        pseudopotential = {"C": "C.UPF"}

        c_scf_params = {
            "input_data": {
                "control": {
                    "calculation": "scf",
                    "restart_mode": "from_scratch",
                    "tprnfor": True,
                    "tstress": True,
                },
                "system": {"ecutwfc": 45, "occupations": "fixed"},
                "electrons": {
                    "diagonalization": "david",
                    "mixing_beta": 0.7,
                    "conv_thr": 1.0e-9,
                },
            },
            "pseudopotentials": pseudopotential,
            "kpts": (3, 3, 3),
        }
        c_scf_results = relax_job(atoms, **c_scf_params)

        c_ph_params = {
            "input_data": {
                "inputph": {
                    "fildyn": "diam.dyn",
                    "fildvscf": "dvscf",
                    "ldisp": True,
                    "nq1": 2,
                    "nq2": 2,
                    "nq3": 2,
                    "tr2_ph": 1.0e-10,
                }
            }
        }

        with caplog.at_level(WARNING):
            c_ph_results = phonon_job(c_scf_results["dir_name"], **c_ph_params)
            assert "Overwriting key 'fildyn'" in caplog.text
            caplog.clear()

        q2r_params = {
            "input_data": {
                "input": {"fildyn": "diam.dyn", "zasr": "crystal", "flfrc": "diam.ifc"}
            }
        }
        q2r_results = q2r_job(c_ph_results["dir_name"], **q2r_params)

        assert q2r_results["parameters"]["input_data"]["input"]["flfrc"] == "q2r.fc"
        assert q2r_results["parameters"]["input_data"]["input"]["fildyn"] == "matdyn"

        dvscf_q2r_params = {
            "input_data": {
                "input": {
                    "prefix": "diam",
                    "fildyn": "diam.dyn",
                    "fildvscf": "dvscf",
                    "wpot_dir": "wpot/",
                    "do_long_range": False,
                    "do_charge_neutral": False,
                }
            }
        }

        with caplog.at_level(WARNING):
            dvscf_q2r_results = dvscf_q2r_job(
                c_ph_results["dir_name"], **dvscf_q2r_params
            )

            assert "Overwriting key 'fildyn'" in caplog.text
            assert "Overwriting key 'wpot_dir'" in caplog.text
            assert "Overwriting key 'prefix'" in caplog.text
            caplog.clear()

        assert (
            dvscf_q2r_results["parameters"]["input_data"]["input"]["fildyn"] == "matdyn"
        )
        assert (
            dvscf_q2r_results["parameters"]["input_data"]["input"]["prefix"] == "pwscf"
        )

        c_nscf_params = {
            "input_data": {
                "control": {"calculation": "nscf", "restart_mode": "from_scratch"},
                "system": {
                    "ecutwfc": 45,
                    "occupations": "fixed",
                    "nbnd": 15,
                    "nosym": True,
                    "noinv": True,
                },
                "electrons": {
                    "diago_full_acc": True,
                    "diagonalization": "david",
                    "mixing_beta": 0.7,
                    "conv_thr": 1.0e-9,
                },
            },
            "pseudopotentials": pseudopotential,
            "kpts": {"path": [[0.0, 0.0, 0.0], [0.365, 0.365, 0.0]]},
        }

        c_nscf_results = non_scf_job(
            c_scf_results["atoms"], c_scf_results["dir_name"], **c_nscf_params
        )

        c_ahc_coarse_params = {
            "input_data": {
                "inputph": {
                    "prefix": "diam",
                    "fildyn": "dyn_dir_ahc_coarse/diam.dyn",
                    "tr2_ph": 1.0e-10,
                    "lDisp": True,
                    "nq1": 3,
                    "nQ2": 3,
                    "nq3": 3,
                    "ldvscf_interpolate": True,
                    "wpOt_dir": "wpot/",
                    "trans": False,
                    "electron_phonon": "ahc",
                    "ahc_nbnd": 8,
                    "AHC_dir": "ahc_dir_coarse",
                }
            }
        }
        c_ahc_coarse_results = phonon_job(
            [dvscf_q2r_results["dir_name"], c_nscf_results["dir_name"]],
            **c_ahc_coarse_params,
        )

        matdyn_coarse_params = {
            "input_data": {
                "input": {
                    "asr": "crystal",
                    "amaSS(1)": 12.01078,
                    "flfrc": "diam.ifc",
                    "flvec": "diam.modes_coarse",
                    "flfrq": "",
                    "fleig": "",
                    "Q_in_band_form": False,
                    "Q_in_cryst_coord": False,
                    "loto_disable": True,
                }
            },
            "additional_cards": [
                "27",
                "0.000000000000000E+00 0.000000000000000E+00 0.000000000000000E+00",
                "-0.333333333333333E+00 0.333333333333333E+00 -0.333333333333333E+00",
                "0.333333333333333E+00 -0.333333333333333E+00 0.333333333333333E+00",
                "0.333333333333333E+00 0.333333333333333E+00 0.333333333333333E+00",
                "0.000000000000000E+00 0.666666666666667E+00 0.000000000000000E+00",
                "0.666666666666667E+00 -0.555111512312578E-16 0.666666666666667E+00",
                "-0.333333333333333E+00 -0.333333333333333E+00 -0.333333333333333E+00",
                "-0.666666666666667E+00 -0.555111512312578E-16 -0.666666666666667E+00",
                "0.000000000000000E+00 -0.666666666666667E+00 0.000000000000000E+00",
                "-0.333333333333333E+00 -0.333333333333333E+00 0.333333333333333E+00",
                "-0.666666666666667E+00 0.000000000000000E+00 0.000000000000000E+00",
                "0.555111512312578E-16 -0.666666666666667E+00 0.666666666666667E+00",
                "0.000000000000000E+00 0.000000000000000E+00 0.666666666666667E+00",
                "-0.333333333333333E+00 0.333333333333333E+00 0.333333333333333E+00",
                "0.333333333333333E+00 -0.333333333333333E+00 0.100000000000000E+01",
                "-0.666666666666667E+00 -0.666666666666667E+00 -0.555111512312578E-16",
                "-0.100000000000000E+01 -0.333333333333333E+00 -0.333333333333333E+00",
                "-0.333333333333333E+00 -0.100000000000000E+01 0.333333333333333E+00",
                "0.333333333333333E+00 0.333333333333333E+00 -0.333333333333333E+00",
                "0.555111512312578E-16 0.666666666666667E+00 -0.666666666666667E+00",
                "0.666666666666667E+00 0.000000000000000E+00 0.000000000000000E+00",
                "0.666666666666667E+00 0.666666666666667E+00 -0.555111512312578E-16",
                "0.333333333333333E+00 0.100000000000000E+01 -0.333333333333333E+00",
                "0.100000000000000E+01 0.333333333333333E+00 0.333333333333333E+00",
                "0.000000000000000E+00 0.000000000000000E+00 -0.666666666666667E+00",
                "-0.333333333333333E+00 0.333333333333333E+00 -0.100000000000000E+01",
                "0.333333333333333E+00 -0.333333333333333E+00 -0.333333333333333E+00",
            ],
        }
        matdyn_coarse_results = matdyn_job(
            [q2r_results["dir_name"], c_ahc_coarse_results["dir_name"]],
            **matdyn_coarse_params,
        )

        postahc_coarse_params = {
            "input_data": {
                "input": {
                    "ahc_dir": "ahc_dir_coarse/",
                    "flvec": "diam.modes_coarse",
                    "nk": 2,
                    "nbnd": 15,
                    "ahc_nbnd": 8,
                    "ahc_nbndskip": 0,
                    "NAT": 2,
                    "NQ": 27,
                    "eta": 0.01,
                    "efermi": 0.4766,
                    "temp_kelvin": 300.0,
                    "amass_amu(1)": 12.01078,
                    "amass_amu(2)": 12.01078,
                }
            }
        }
        postahc_coarse_results = postahc_job(
            [c_ahc_coarse_results["dir_name"], matdyn_coarse_results["dir_name"]],
            **postahc_coarse_params,
        )

        assert Path(postahc_coarse_results["dir_name"], "postahc.out.gz").exists()
        assert Path(postahc_coarse_results["dir_name"], "ahc_dir").exists()
        assert Path(postahc_coarse_results["dir_name"], "selfen_real.dat.gz").exists()


def test_phonon_dvscf_q2r_inplace(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    with change_settings({"ESPRESSO_PSEUDO": tmp_path}):
        copy_decompress_files(DATA_DIR, ["C.UPF.gz"], tmp_path)

        atoms = bulk("C")

        pseudopotentials = {"C": "C.UPF"}

        c_scf_results = static_job(
            atoms, pseudopotentials=pseudopotentials, kspacing=0.1
        )

        c_ph_params = {
            "input_data": {
                "inputph": {
                    "fildvscf": "dvscf",
                    "ldisp": True,
                    "nq1": 2,
                    "nq2": 2,
                    "nq3": 2,
                    "tr2_ph": 1.0e-6,
                }
            }
        }

        c_ph_results = phonon_job(prev_outdir=c_scf_results["dir_name"], **c_ph_params)

        dvscf_q2r_params = {
            "input_data": {
                "input": {
                    "fildvscf": "dvscf",
                    "do_long_range": False,
                    "do_charge_neutral": False,
                }
            }
        }

        dvscf_q2r_results = dvscf_q2r_job(
            prev_outdir=c_scf_results["dir_name"],
            copy_files=Copy(src_dir=c_ph_results["dir_name"]),
            **dvscf_q2r_params,
        )

        assert (
            dvscf_q2r_results["parameters"]["input_data"]["input"]["outdir"]
            == c_scf_results["dir_name"]
        )
        assert Path(dvscf_q2r_results["dir_name"], "matdyn0.gz").exists()
