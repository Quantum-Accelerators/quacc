from __future__ import annotations

from shutil import which

import pytest

from quacc import SETTINGS

pytestmark = pytest.mark.skipif(
    which(str(SETTINGS.ESPRESSO_BIN_DIR / SETTINGS.ESPRESSO_BINARIES["pw"])) is None
    or which(str(SETTINGS.ESPRESSO_BIN_DIR / SETTINGS.ESPRESSO_BINARIES["ph"])) is None,
    reason="QE not installed",
)

from pathlib import Path

from ase.build import bulk
from ase.io.espresso import read_fortran_namelist
from monty.shutil import decompress_file
from numpy.testing import assert_allclose, assert_array_equal

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

DEFAULT_SETTINGS = SETTINGS.model_copy()
DATA_DIR = Path(__file__).parent / "data"


def test_phonon_job(tmp_path, monkeypatch, ESPRESSO_PARALLEL_INFO):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    atoms = bulk("Li")

    copy_decompress_files(DATA_DIR, ["Li.upf.gz"], tmp_path)

    SETTINGS.ESPRESSO_PSEUDO = tmp_path

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
        parallel_info=ESPRESSO_PARALLEL_INFO,
    )

    ph_results = phonon_job(
        pw_results["dir_name"],
        input_data=ph_loose,
        parallel_info=ESPRESSO_PARALLEL_INFO,
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

    SETTINGS.ESPRESSO_PSEUDO = DEFAULT_SETTINGS.ESPRESSO_PSEUDO


def test_phonon_job_list_to_do(tmp_path, monkeypatch, ESPRESSO_PARALLEL_INFO):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    atoms = bulk("Li")

    copy_decompress_files(DATA_DIR, ["Li.upf.gz"], tmp_path)

    SETTINGS.ESPRESSO_PSEUDO = tmp_path

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
        parallel_info=ESPRESSO_PARALLEL_INFO,
    )

    qpts = [(0, 0, 0, 1), (1 / 3, 0, 0, 1), (1 / 2, 0, 0, 1)]

    nat_todo = [1]

    ph_results = phonon_job(
        pw_results["dir_name"],
        input_data=ph_loose,
        qpts=qpts,
        nat_todo_indices=nat_todo,
        parallel_info=ESPRESSO_PARALLEL_INFO,
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

    SETTINGS.ESPRESSO_PSEUDO = DEFAULT_SETTINGS.ESPRESSO_PSEUDO


def test_q2r_job(tmp_path, monkeypatch, ESPRESSO_PARALLEL_INFO):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    copy_decompress_files(DATA_DIR / "q2r_test", "matdyn", tmp_path)

    additional_cards = ["1 1 1", "1", "matdyn"]

    q2r_results = q2r_job(
        tmp_path,
        additional_cards=additional_cards,
        parallel_info=ESPRESSO_PARALLEL_INFO,
    )

    assert Path(q2r_results["dir_name"], "q2r.fc.gz").exists()

    decompress_file(Path(q2r_results["dir_name"], "q2r.in.gz"))

    with Path(q2r_results["dir_name"], "q2r.in").open() as f:
        recycled_input = read_fortran_namelist(f)

    assert recycled_input[0]["input"].pop("flfrc") == "q2r.fc"

    assert recycled_input[0]["input"] == {"fildyn": "matdyn"}
    assert recycled_input[1] == [*additional_cards, "EOF"]


def test_matdyn_job(tmp_path, monkeypatch, ESPRESSO_PARALLEL_INFO):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    copy_decompress_files(DATA_DIR / "matdyn_test", "q2r.fc", tmp_path)

    input_data = {"input": {"dos": True, "nk1": 4, "nk2": 4, "nk3": 4}}
    matdyn_results = matdyn_job(
        tmp_path, input_data=input_data, parallel_info=ESPRESSO_PARALLEL_INFO
    )

    assert Path(matdyn_results["dir_name"], "q2r.fc.gz").exists()
    assert Path(matdyn_results["dir_name"], "matdyn.dos.gz").exists()
    assert Path(matdyn_results["dir_name"], "matdyn.freq.gz").exists()
    assert Path(matdyn_results["dir_name"], "matdyn.modes.gz").exists()
    assert matdyn_results["results"]["matdyn_results"]["phonon_dos"].shape == (561, 3)

    decompress_file(Path(matdyn_results["dir_name"], "matdyn.in.gz"))

    with Path(matdyn_results["dir_name"], "matdyn.in").open() as f:
        recycled_input = read_fortran_namelist(f)

    assert recycled_input[0]["input"].pop("flfrc") == "q2r.fc"

    assert recycled_input[0]["input"]["fldos"] == "matdyn.dos"
    assert recycled_input[0]["input"]["flfrq"] == "matdyn.freq"
    assert recycled_input[0]["input"]["flvec"] == "matdyn.modes"
    assert recycled_input[0]["input"]["fleig"] == "matdyn.eig"


def test_phonon_dos_flow(tmp_path, monkeypatch, ESPRESSO_PARALLEL_INFO):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    atoms = bulk("Li")

    copy_decompress_files(DATA_DIR, ["Li.upf.gz"], tmp_path)

    SETTINGS.ESPRESSO_PSEUDO = tmp_path

    input_data = {
        "control": {"calculation": "scf"},
        "system": {"occupations": "smearing", "smearing": "cold", "degauss": 0.02},
        "electrons": {"mixing_mode": "TF", "mixing_beta": 0.7, "conv_thr": 1.0e-6},
    }

    pseudopotentials = {"Li": "Li.upf"}

    job_params = {
        "relax_job": {
            "pseudopotentials": pseudopotentials,
            "input_data": input_data,
            "kspacing": 1.0,
            "parallel_info": ESPRESSO_PARALLEL_INFO,
        },
        "phonon_job": {"parallel_info": ESPRESSO_PARALLEL_INFO},
        "q2r_job": {"parallel_info": ESPRESSO_PARALLEL_INFO},
        "matdyn_job": {"parallel_info": ESPRESSO_PARALLEL_INFO},
    }

    assert phonon_dos_flow(atoms, job_params=job_params)

    SETTINGS.ESPRESSO_PSEUDO = DEFAULT_SETTINGS.ESPRESSO_PSEUDO


def test_phonon_calculation_spin_orbit_example_06(
    tmp_path, monkeypatch, ESPRESSO_PARALLEL_INFO
):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    # Create input structure
    pt_atoms = bulk("Pt")

    copy_decompress_files(DATA_DIR, ["Pt.rel-pz-n-rrkjus.UPF.gz"], tmp_path)

    SETTINGS.ESPRESSO_PSEUDO = tmp_path
    # Set up pseudopotential
    pseudopotential = {"Pt": "Pt.rel-pz-n-rrkjus.UPF"}

    # Pt calculations
    pt_relax_params = {
        "input_data": {
            "control": {
                "calculation": "scf",
                "restart_mode": "from_scratch",
                "tprnfor": True,
                "tstress": True,
            },
            "system": {
                "lspinorb": True,
                "noncolin": True,
                "starting_magnetization": 0.0,
                "occupations": "smearing",
                "degauss": 0.02,
                "smearing": "mv",
                "ecutwfc": 30.0,
                "ecutrho": 250.0,
            },
            "electrons": {"mixing_beta": 0.7, "conv_thr": 1.0e-8},
        },
        "pseudopotentials": pseudopotential,
        "kpts": (2, 2, 2),
        "parallel_info": ESPRESSO_PARALLEL_INFO,
    }
    pt_relax_results = relax_job(pt_atoms, **pt_relax_params)

    pt_phonon_params = {
        "input_data": {
            "inputph": {"amass(1)": 195.078, "fildyn": "ptdyn", "tr2_ph": 1.0e-16}
        },
        "parallel_info": ESPRESSO_PARALLEL_INFO,
    }
    pt_phonon_results = phonon_job(
        pt_relax_results["dir_name"], **pt_phonon_params, qpts=(0.0, 0.0, 0.0)
    )

    pt_phonon_x_params = {
        "input_data": {
            "inputph": {"amass(1)": 195.078, "fildyn": "ptdyn", "tr2_ph": 1.0e-16}
        },
        "parallel_info": ESPRESSO_PARALLEL_INFO,
    }
    pt_phonon_x_results = phonon_job(
        pt_relax_results["dir_name"], **pt_phonon_x_params, qpts=(1.0, 0.0, 0.0)
    )

    pt_phonon_dispersions_params = {
        "input_data": {
            "inputph": {
                "amass(1)": 195.078,
                "fildyn": "ptdyn",
                "tr2_ph": 1.0e-16,
                "ldisp": True,
                "nq1": 4,
                "nq2": 4,
                "nq3": 4,
            }
        },
        "parallel_info": ESPRESSO_PARALLEL_INFO,
    }
    pt_phonon_dispersions_results = phonon_job(
        pt_relax_results["dir_name"], **pt_phonon_dispersions_params
    )

    SETTINGS.ESPRESSO_PSEUDO = DEFAULT_SETTINGS.ESPRESSO_PSEUDO


def test_phonon_calculation_si_spin_orbit(
    tmp_path, monkeypatch, ESPRESSO_PARALLEL_INFO
):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    SETTINGS.ESPRESSO_PSEUDO = tmp_path

    copy_decompress_files(DATA_DIR, ["Si.rel-pbe.rrkj.UPF.gz"], tmp_path)

    si_atoms = bulk("Si")

    pseudopotential = {"Si": "Si.rel-pbe.rrkj.UPF"}

    si_relax_params = {
        "input_data": {
            "control": {"calculation": "scf", "restart_mode": "from_scratch"},
            "system": {"ecutwfc": 20.0, "noncolin": True, "lspinorb": True},
            "electrons": {"mixing_beta": 0.7, "conv_thr": 1.0e-10},
        },
        "pseudopotentials": pseudopotential,
        "kpts": (2, 2, 2),
        "parallel_info": ESPRESSO_PARALLEL_INFO,
    }
    si_relax_results = relax_job(si_atoms, **si_relax_params)

    si_phonon_params = {
        "input_data": {
            "inputph": {
                "tr2_ph": 1.0e-16,
                "epsil": True,
                "fildyn": "Sig.dyn",
                "amass(1)": 28.0855,
            }
        },
        "parallel_info": ESPRESSO_PARALLEL_INFO,
    }
    si_phonon_results = phonon_job(
        si_relax_results["dir_name"], **si_phonon_params, qpts=(0.0, 0.0, 0.0)
    )

    assert si_relax_results is not None
    assert si_phonon_results is not None

    SETTINGS.ESPRESSO_PSEUDO = DEFAULT_SETTINGS.ESPRESSO_PSEUDO


def test_phonon_calculation_c_noncollinear(
    tmp_path, monkeypatch, ESPRESSO_PARALLEL_INFO
):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    SETTINGS.ESPRESSO_PSEUDO = tmp_path

    copy_decompress_files(DATA_DIR, ["C.pz-rrkjus.UPF.gz"], tmp_path)

    c_atoms = bulk("C")

    pseudopotential = {"C": "C.pz-rrkjus.UPF"}

    c_relax_params = {
        "input_data": {
            "control": {"calculation": "scf", "restart_mode": "from_scratch"},
            "system": {"ecutwfc": 27.0, "ecutrho": 300.0, "noncolin": True},
            "electrons": {"mixing_beta": 0.7, "conv_thr": 1.0e-9},
        },
        "pseudopotentials": pseudopotential,
        "kpts": (4, 4, 4),
        "parallel_info": ESPRESSO_PARALLEL_INFO,
    }
    c_relax_results = relax_job(c_atoms, **c_relax_params)

    c_phonon_params = {
        "input_data": {"inputph": {"tr2_ph": 1.0e-14, "epsil": True}},
        "parallel_info": ESPRESSO_PARALLEL_INFO,
    }
    c_phonon_results = phonon_job(
        c_relax_results["dir_name"], **c_phonon_params, qpts=(0.0, 0.0, 0.0)
    )

    assert c_relax_results is not None
    assert c_phonon_results is not None

    SETTINGS.ESPRESSO_PSEUDO = DEFAULT_SETTINGS.ESPRESSO_PSEUDO


def test_phonon_induced_renormalization(tmp_path, monkeypatch, ESPRESSO_PARALLEL_INFO):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    SETTINGS.ESPRESSO_PSEUDO = tmp_path

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
            "system": {"ecutwfc": 60, "occupations": "fixed"},
            "electrons": {
                "diagonalization": "david",
                "mixing_beta": 0.7,
                "conv_thr": 1.0e-12,
            },
        },
        "pseudopotentials": pseudopotential,
        "kpts": (6, 6, 6),
    }
    c_scf_results = relax_job(
        atoms, **c_scf_params, parallel_info=ESPRESSO_PARALLEL_INFO
    )

    c_ph_params = {
        "input_data": {
            "inputph": {
                "fildyn": "diam.dyn",
                "fildvscf": "dvscf",
                "ldisp": True,
                "nq1": 3,
                "nq2": 3,
                "nq3": 3,
                "tr2_ph": 1.0e-16,
            }
        }
    }
    c_ph_results = phonon_job(
        c_scf_results["dir_name"], **c_ph_params, parallel_info=ESPRESSO_PARALLEL_INFO
    )

    q2r_params = {
        "input_data": {
            "input": {"fildyn": "diam.dyn", "zasr": "crystal", "flfrc": "diam.ifc"}
        }
    }
    q2r_results = q2r_job(
        c_ph_results["dir_name"], **q2r_params, parallel_info=ESPRESSO_PARALLEL_INFO
    )

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
        },
        "parallel_info": ESPRESSO_PARALLEL_INFO,
    }
    dvscf_q2r_results = dvscf_q2r_job(
        c_ph_results["dir_name"],
        **dvscf_q2r_params,
        parallel_info=ESPRESSO_PARALLEL_INFO,
    )

    c_nscf_params = {
        "input_data": {
            "control": {"calculation": "nscf", "restart_mode": "from_scratch"},
            "system": {
                "ecutwfc": 60,
                "occupations": "fixed",
                "nbnd": 15,
                "nosym": True,
                "noinv": True,
            },
            "electrons": {
                "diago_full_acc": True,
                "diagonalization": "david",
                "mixing_beta": 0.7,
                "conv_thr": 1.0e-10,
            },
        },
        "pseudopotentials": pseudopotential,
        "kpts": {"path": [[0.0, 0.0, 0.0], [0.365, 0.365, 0.0]]},
        "parallel_info": ESPRESSO_PARALLEL_INFO,
    }

    c_nscf_results = non_scf_job(**c_nscf_params, parallel_info=ESPRESSO_PARALLEL_INFO)

    c_ahc_coarse_params = {
        "input_data": {
            "inputph": {
                "prefix": "diam",
                "fildyn": "dyn_dir_ahc_coarse/diam.dyn",
                "tr2_ph": 1.0e-20,
                "ldisp": True,
                "nq1": 3,
                "nq2": 3,
                "nq3": 3,
                "ldvscf_interpolate": True,
                "wpot_dir": "wpot/",
                "trans": False,
                "electron_phonon": "ahc",
                "ahc_nbnd": 8,
                "ahc_dir": "ahc_dir_coarse",
            }
        },
        "parallel_info": ESPRESSO_PARALLEL_INFO,
    }
    c_ahc_coarse_results = phonon_job(c_nscf_results["dir_name"], **c_ahc_coarse_params)

    matdyn_coarse_params = {
        "input_data": {
            "input": {
                "asr": "crystal",
                "amass(1)": 12.01078,
                "flfrc": "diam.ifc",
                "flvec": "diam.modes_coarse",
                "flfrq": "",
                "fleig": "",
                "q_in_band_form": False,
                "q_in_cryst_coord": False,
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
        "parallel_info": ESPRESSO_PARALLEL_INFO,
    }
    matdyn_coarse_results = matdyn_job(
        c_ahc_coarse_results["dir_name"],
        **matdyn_coarse_params,
        parallel_info=ESPRESSO_PARALLEL_INFO,
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
                "nat": 2,
                "nq": 27,
                "eta": 0.01,
                "efermi": 0.4766,
                "temp_kelvin": 300.0,
                "amass_amu(1)": 12.01078,
                "amass_amu(2)": 12.01078,
            }
        },
        "parallel_info": ESPRESSO_PARALLEL_INFO,
    }
    postahc_coarse_results = postahc_job(
        matdyn_coarse_results["dir_name"],
        **postahc_coarse_params,
        parallel_info=ESPRESSO_PARALLEL_INFO,
    )

    c_ahc_fine_params = {
        "input_data": {
            "inputph": {
                "prefix": "diam",
                "fildyn": "dyn_dir_ahc_fine/diam.dyn",
                "tr2_ph": 1.0e-20,
                "ldisp": True,
                "nq1": 4,
                "nq2": 4,
                "nq3": 4,
                "ldvscf_interpolate": True,
                "wpot_dir": "wpot/",
                "trans": False,
                "electron_phonon": "ahc",
                "ahc_nbnd": 8,
                "ahc_dir": "ahc_dir_fine",
                "skip_upperfan": True,
            }
        },
        "parallel_info": ESPRESSO_PARALLEL_INFO,
    }
    c_ahc_fine_results = phonon_job(
        c_nscf_results["dir_name"],
        **c_ahc_fine_params,
        parallel_info=ESPRESSO_PARALLEL_INFO,
    )

    matdyn_fine_params = {
        "input_data": {
            "input": {
                "asr": "crystal",
                "amass(1)": 12.01078,
                "flfrc": "diam.ifc",
                "flvec": "diam.modes_fine",
                "flfrq": "",
                "fleig": "",
                "q_in_band_form": False,
                "q_in_cryst_coord": False,
                "loto_disable": True,
            }
        },
        "additional_cards": [
            "64",
            "0.000000000000000E+00   0.000000000000000E+00   0.000000000000000E+00",
            "-0.250000000000000E+00   0.250000000000000E+00  -0.250000000000000E+00",
            "0.500000000000000E+00  -0.500000000000000E+00   0.500000000000000E+00",
            "0.250000000000000E+00  -0.250000000000000E+00   0.250000000000000E+00",
            "0.250000000000000E+00   0.250000000000000E+00   0.250000000000000E+00",
            "0.000000000000000E+00   0.500000000000000E+00   0.000000000000000E+00",
            "0.750000000000000E+00  -0.250000000000000E+00   0.750000000000000E+00",
            "0.500000000000000E+00   0.000000000000000E+00   0.500000000000000E+00",
            "-0.500000000000000E+00  -0.500000000000000E+00  -0.500000000000000E+00",
            "-0.750000000000000E+00  -0.250000000000000E+00  -0.750000000000000E+00",
            "0.000000000000000E+00  -0.100000000000000E+01   0.000000000000000E+00",
            "-0.250000000000000E+00  -0.750000000000000E+00  -0.250000000000000E+00",
            "-0.250000000000000E+00  -0.250000000000000E+00  -0.250000000000000E+00",
            "-0.500000000000000E+00   0.000000000000000E+00  -0.500000000000000E+00",
            "0.250000000000000E+00  -0.750000000000000E+00   0.250000000000000E+00",
            "0.000000000000000E+00  -0.500000000000000E+00   0.000000000000000E+00",
            "-0.250000000000000E+00  -0.250000000000000E+00   0.250000000000000E+00",
            "-0.500000000000000E+00   0.000000000000000E+00   0.000000000000000E+00",
            "0.250000000000000E+00  -0.750000000000000E+00   0.750000000000000E+00",
            "0.000000000000000E+00  -0.500000000000000E+00   0.500000000000000E+00",
            "0.000000000000000E+00   0.000000000000000E+00   0.500000000000000E+00",
            "-0.250000000000000E+00   0.250000000000000E+00   0.250000000000000E+00",
            "0.500000000000000E+00  -0.500000000000000E+00   0.100000000000000E+01",
            "0.250000000000000E+00  -0.250000000000000E+00   0.750000000000000E+00",
            "-0.750000000000000E+00  -0.750000000000000E+00  -0.250000000000000E+00",
            "-0.100000000000000E+01  -0.500000000000000E+00  -0.500000000000000E+00",
            "-0.250000000000000E+00  -0.125000000000000E+01   0.250000000000000E+00",
            "-0.500000000000000E+00  -0.100000000000000E+01   0.000000000000000E+00",
            "-0.500000000000000E+00  -0.500000000000000E+00   0.000000000000000E+00",
            "-0.750000000000000E+00  -0.250000000000000E+00  -0.250000000000000E+00",
            "0.000000000000000E+00  -0.100000000000000E+01   0.500000000000000E+00",
            "-0.250000000000000E+00  -0.750000000000000E+00   0.250000000000000E+00",
            "0.500000000000000E+00   0.500000000000000E+00  -0.500000000000000E+00",
            "0.250000000000000E+00   0.750000000000000E+00  -0.750000000000000E+00",
            "0.100000000000000E+01   0.000000000000000E+00   0.000000000000000E+00",
            "0.750000000000000E+00   0.250000000000000E+00  -0.250000000000000E+00",
            "0.750000000000000E+00   0.750000000000000E+00  -0.250000000000000E+00",
            "0.500000000000000E+00   0.100000000000000E+01  -0.500000000000000E+00",
            "0.125000000000000E+01   0.250000000000000E+00   0.250000000000000E+00",
            "0.100000000000000E+01   0.500000000000000E+00   0.000000000000000E+00",
            "0.000000000000000E+00   0.000000000000000E+00  -0.100000000000000E+01",
            "-0.250000000000000E+00   0.250000000000000E+00  -0.125000000000000E+01",
            "0.500000000000000E+00  -0.500000000000000E+00  -0.500000000000000E+00",
            "0.250000000000000E+00  -0.250000000000000E+00  -0.750000000000000E+00",
            "0.250000000000000E+00   0.250000000000000E+00  -0.750000000000000E+00",
            "0.000000000000000E+00   0.500000000000000E+00  -0.100000000000000E+01",
            "0.750000000000000E+00  -0.250000000000000E+00  -0.250000000000000E+00",
            "0.500000000000000E+00   0.000000000000000E+00  -0.500000000000000E+00",
            "0.250000000000000E+00   0.250000000000000E+00  -0.250000000000000E+00",
            "0.000000000000000E+00   0.500000000000000E+00  -0.500000000000000E+00",
            "0.750000000000000E+00  -0.250000000000000E+00   0.250000000000000E+00",
            "0.500000000000000E+00   0.000000000000000E+00   0.000000000000000E+00",
            "0.500000000000000E+00   0.500000000000000E+00   0.000000000000000E+00",
            "0.250000000000000E+00   0.750000000000000E+00  -0.250000000000000E+00",
            "0.100000000000000E+01   0.000000000000000E+00   0.500000000000000E+00",
            "0.750000000000000E+00   0.250000000000000E+00   0.250000000000000E+00",
            "-0.250000000000000E+00  -0.250000000000000E+00  -0.750000000000000E+00",
            "-0.500000000000000E+00   0.000000000000000E+00  -0.100000000000000E+01",
            "0.250000000000000E+00  -0.750000000000000E+00  -0.250000000000000E+00",
            "0.000000000000000E+00  -0.500000000000000E+00  -0.500000000000000E+00",
            "0.000000000000000E+00   0.000000000000000E+00  -0.500000000000000E+00",
            "-0.250000000000000E+00   0.250000000000000E+00  -0.750000000000000E+00",
            "0.500000000000000E+00  -0.500000000000000E+00   0.000000000000000E+00",
            "0.250000000000000E+00  -0.250000000000000E+00  -0.250000000000000E+00",
        ],
        "parallel_info": ESPRESSO_PARALLEL_INFO,
    }
    matdyn_fine_results = matdyn_job(
        c_ahc_fine_results["dir_name"], **matdyn_fine_params
    )

    postahc_fine_params = {
        "input_data": {
            "input": {
                "ahc_dir": "ahc_dir_fine/",
                "flvec": "diam.modes_fine",
                "nk": 2,
                "nbnd": 15,
                "ahc_nbnd": 8,
                "ahc_nbndskip": 0,
                "nat": 2,
                "nq": 64,
                "eta": 0.01,
                "efermi": 0.4766,
                "temp_kelvin": 300.0,
                "skip_upperfan": True,
                "skip_dw": True,
                "amass_amu(1)": 12.01078,
                "amass_amu(2)": 12.01078,
            }
        },
        "parallel_info": ESPRESSO_PARALLEL_INFO,
    }
    postahc_fine_results = postahc_job(**postahc_fine_params)

    assert c_scf_results is not None
    assert c_ph_results is not None
    assert q2r_results is not None
    assert dvscf_q2r_results is not None
    assert c_nscf_results is not None
    assert c_ahc_coarse_results is not None
    assert matdyn_coarse_results is not None
    assert postahc_coarse_results is not None
    assert c_ahc_fine_results is not None
    assert matdyn_fine_results is not None
    assert postahc_fine_results is not None

    SETTINGS.ESPRESSO_PSEUDO = DEFAULT_SETTINGS.ESPRESSO_PSEUDO
