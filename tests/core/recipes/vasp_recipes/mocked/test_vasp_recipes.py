import pytest
from ase.build import bulk, molecule

from quacc import SETTINGS
from quacc.recipes.vasp.core import double_relax_job, relax_job, static_job
from quacc.recipes.vasp.mp import (
    mp_gga_relax_flow,
    mp_gga_relax_job,
    mp_gga_static_job,
    mp_metagga_prerelax_job,
    mp_metagga_relax_flow,
    mp_metagga_relax_job,
    mp_metagga_static_job,
)
from quacc.recipes.vasp.qmof import qmof_relax_job
from quacc.recipes.vasp.slabs import bulk_to_slabs_flow
from quacc.recipes.vasp.slabs import relax_job as slab_relax_job
from quacc.recipes.vasp.slabs import slab_to_ads_flow
from quacc.recipes.vasp.slabs import static_job as slab_static_job

DEFAULT_SETTINGS = SETTINGS.model_copy()


def test_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Al")

    output = static_job(atoms)
    assert output["nsites"] == len(atoms)
    assert "isym" not in output["parameters"]
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["lwave"] is True
    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["efermi"] == "midgap"

    output = static_job(atoms, ncore=2, kpar=4)
    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["ncore"] == 2
    assert output["parameters"]["kpar"] == 4

    output = static_job(atoms, preset="QMOFSet", ismear=0, sigma=0.01, nedos=None)
    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["ismear"] == 0
    assert output["parameters"]["sigma"] == 0.01

    output = static_job(atoms, ivdw=11, lasph=False, prec=None, lwave=None, efermi=None)
    assert output["parameters"]["ivdw"] == 11
    assert output["parameters"]["lasph"] is False
    assert "prec" not in output["parameters"]
    assert "lwave" not in output["parameters"]
    assert "efermi" not in output["parameters"]


def test_static_job_incar_copilot_aggressive(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Al")

    SETTINGS.VASP_INCAR_COPILOT = "aggressive"
    output = static_job(atoms, ivdw=11, lasph=False, prec=None, lwave=None, efermi=None)
    assert output["parameters"]["ivdw"] == 11
    assert output["parameters"]["lasph"] is False
    assert "prec" not in output["parameters"]
    assert "lwave" not in output["parameters"]
    assert output["parameters"]["efermi"] == "midgap"
    SETTINGS.VASP_INCAR_COPILOT = DEFAULT_SETTINGS.VASP_INCAR_COPILOT


def test_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Al")

    output = relax_job(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isif"] == 3
    assert output["parameters"]["lwave"] is False
    assert output["parameters"]["encut"] == 520

    output = relax_job(atoms, nelmin=6)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isif"] == 3
    assert output["parameters"]["lwave"] is False
    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["nelmin"] == 6

    output = relax_job(atoms, relax_cell=False)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["lwave"] is False
    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["isif"] == 2


def test_doublerelax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Al")

    output = double_relax_job(atoms)
    assert output["relax1"]["nsites"] == len(atoms)
    assert output["relax1"]["parameters"]["isym"] == 0
    assert output["relax1"]["parameters"]["nsw"] > 0
    assert output["relax1"]["parameters"]["isif"] == 3
    assert output["relax1"]["parameters"]["lwave"] is False
    assert output["relax1"]["parameters"]["encut"] == 520
    assert output["relax2"]["nsites"] == len(atoms)
    assert output["relax2"]["parameters"]["isym"] == 0
    assert output["relax2"]["parameters"]["nsw"] > 0
    assert output["relax2"]["parameters"]["isif"] == 3
    assert output["relax2"]["parameters"]["lwave"] is False
    assert output["relax2"]["parameters"]["encut"] == 520

    output = double_relax_job(atoms, relax2_kwargs={"nelmin": 6})
    assert output["relax1"]["nsites"] == len(atoms)
    assert output["relax1"]["parameters"]["isym"] == 0
    assert output["relax1"]["parameters"]["nsw"] > 0
    assert output["relax1"]["parameters"]["isif"] == 3
    assert output["relax1"]["parameters"]["lwave"] is False
    assert output["relax1"]["parameters"]["encut"] == 520
    assert output["relax2"]["nsites"] == len(atoms)
    assert output["relax2"]["parameters"]["isym"] == 0
    assert output["relax2"]["parameters"]["nsw"] > 0
    assert output["relax2"]["parameters"]["isif"] == 3
    assert output["relax2"]["parameters"]["lwave"] is False
    assert output["relax2"]["parameters"]["encut"] == 520
    assert output["relax2"]["parameters"]["nelmin"] == 6

    output = double_relax_job(atoms, relax_cell=False)
    assert output["relax1"]["nsites"] == len(atoms)
    assert output["relax1"]["parameters"]["isym"] == 0
    assert output["relax1"]["parameters"]["nsw"] > 0
    assert output["relax1"]["parameters"]["lwave"] is False
    assert output["relax1"]["parameters"]["isif"] == 2
    assert output["relax1"]["parameters"]["encut"] == 520
    assert output["relax2"]["nsites"] == len(atoms)
    assert output["relax2"]["parameters"]["isym"] == 0
    assert output["relax2"]["parameters"]["nsw"] > 0
    assert output["relax2"]["parameters"]["lwave"] is False
    assert output["relax2"]["parameters"]["encut"] == 520
    assert output["relax2"]["parameters"]["isif"] == 2

    assert double_relax_job(atoms, relax1_kwargs={"kpts": [1, 1, 1]})


def test_slab_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Al")

    output = slab_static_job(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["idipol"] == 3
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["lvhar"] is True
    assert output["parameters"]["encut"] == 450

    output = slab_static_job(atoms, nelmin=6)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["idipol"] == 3
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["lvhar"] is True
    assert output["parameters"]["encut"] == 450
    assert output["parameters"]["nelmin"] == 6

    output = slab_static_job(atoms, encut=None)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["idipol"] == 3
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["lvhar"] is True
    assert "encut" not in output["parameters"]


def test_slab_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Al")

    output = slab_relax_job(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isif"] == 2
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["lwave"] is False
    assert output["parameters"]["encut"] == 450

    output = slab_relax_job(atoms, nelmin=6)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isif"] == 2
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["lwave"] is False
    assert output["parameters"]["encut"] == 450
    assert output["parameters"]["nelmin"] == 6


def test_slab_dynamic_jobs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Al")

    ### --------- Test bulk_to_slabs_flow --------- ###

    outputs = bulk_to_slabs_flow(atoms, run_static=False)
    assert len(outputs) == 4
    assert outputs[0]["nsites"] == 45
    assert outputs[1]["nsites"] == 45
    assert outputs[2]["nsites"] == 54
    assert outputs[3]["nsites"] == 42
    assert [output["parameters"]["isif"] == 2 for output in outputs]

    outputs = bulk_to_slabs_flow(atoms)
    assert len(outputs) == 4
    assert outputs[0]["nsites"] == 45
    assert outputs[1]["nsites"] == 45
    assert outputs[2]["nsites"] == 54
    assert outputs[3]["nsites"] == 42
    assert [output["parameters"]["nsw"] == 0 for output in outputs]

    outputs = bulk_to_slabs_flow(
        atoms,
        job_params={"relax_job": {"preset": "SlabSet", "nelmin": 6}},
        run_static=False,
    )
    assert len(outputs) == 4
    assert outputs[0]["nsites"] == 45
    assert outputs[1]["nsites"] == 45
    assert outputs[2]["nsites"] == 54
    assert outputs[3]["nsites"] == 42
    assert [output["parameters"]["isif"] == 2 for output in outputs]
    assert [output["parameters"]["nelmin"] == 6 for output in outputs]
    assert [output["parameters"]["encut"] == 450 for output in outputs]

    outputs = bulk_to_slabs_flow(
        atoms, job_params={"relax_job": {"preset": "SlabSet", "nelmin": 6}}
    )
    assert len(outputs) == 4
    assert outputs[0]["nsites"] == 45
    assert outputs[1]["nsites"] == 45
    assert outputs[2]["nsites"] == 54
    assert outputs[3]["nsites"] == 42
    assert [output["parameters"]["nsw"] == 0 for output in outputs]
    assert [output["parameters"]["nelmin"] == 6 for output in outputs]
    assert [output["parameters"]["encut"] == 450 for output in outputs]

    ### --------- Test slab_to_ads_flow --------- ###
    atoms = outputs[0]["atoms"]
    adsorbate = molecule("H2")

    outputs = slab_to_ads_flow(atoms, adsorbate, run_static=False)

    assert [output["nsites"] == 82 for output in outputs]
    assert [output["parameters"]["isif"] == 2 for output in outputs]

    outputs = slab_to_ads_flow(atoms, adsorbate)
    assert [output["nsites"] == 82 for output in outputs]
    assert [output["parameters"]["nsw"] == 0 for output in outputs]

    outputs = slab_to_ads_flow(
        atoms,
        adsorbate,
        job_params={"relax_job": {"preset": "SlabSet", "nelmin": 6}},
        run_static=False,
    )

    assert [output["nsites"] == 82 for output in outputs]
    assert [output["parameters"]["isif"] == 2 for output in outputs]
    assert [output["parameters"]["nelmin"] == 6 for output in outputs]
    assert [output["parameters"]["encut"] == 450 for output in outputs]

    outputs = slab_to_ads_flow(
        atoms, adsorbate, job_params={"relax_job": {"preset": "SlabSet", "nelmin": 6}}
    )

    assert [output["nsites"] == 82 for output in outputs]
    assert [output["parameters"]["nsw"] == 0 for output in outputs]
    assert [output["parameters"]["nelmin"] == 6 for output in outputs]
    assert [output["parameters"]["encut"] == 450 for output in outputs]

    adsorbate2 = molecule("CH3")
    adsorbate2.set_initial_magnetic_moments([1, 0, 0, 0])
    outputs = slab_to_ads_flow(atoms, adsorbate2)
    assert [output["nsites"] == 84 for output in outputs]
    assert [output["parameters"]["nsw"] == 0 for output in outputs]


def test_qmof(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Al")
    output = qmof_relax_job(atoms)
    assert output["prerelax_lowacc"]["nsites"] == len(atoms)
    assert output["prerelax_lowacc"]["parameters"]["sigma"] == 0.01
    assert output["prerelax_lowacc"]["parameters"]["isym"] == 0
    assert output["prerelax_lowacc"]["parameters"]["nsw"] == 0
    assert "isif" not in output["prerelax_lowacc"]["parameters"]
    assert "encut" not in output["prerelax_lowacc"]["parameters"]

    assert output["position_relax_lowacc"]["nsites"] == len(atoms)
    assert output["position_relax_lowacc"]["parameters"]["sigma"] == 0.01
    assert output["position_relax_lowacc"]["parameters"]["isym"] == 0
    assert output["position_relax_lowacc"]["parameters"]["nsw"] > 0
    assert output["position_relax_lowacc"]["parameters"]["isif"] == 2
    assert "encut" not in output["prerelax_lowacc"]["parameters"]

    assert output["volume_relax_lowacc"]["nsites"] == len(atoms)
    assert output["volume_relax_lowacc"]["parameters"]["encut"] == 520
    assert output["volume_relax_lowacc"]["parameters"]["sigma"] == 0.01
    assert output["volume_relax_lowacc"]["parameters"]["isym"] == 0
    assert output["volume_relax_lowacc"]["parameters"]["nsw"] > 0
    assert output["volume_relax_lowacc"]["parameters"]["isif"] == 3

    assert output["double_relax"][0]["nsites"] == len(atoms)
    assert output["double_relax"][0]["parameters"]["encut"] == 520
    assert output["double_relax"][0]["parameters"]["sigma"] == 0.01
    assert output["double_relax"][0]["parameters"]["isym"] == 0
    assert output["double_relax"][0]["parameters"]["nsw"] > 0
    assert output["double_relax"][0]["parameters"]["isif"] == 3

    assert output["double_relax"][1]["nsites"] == len(atoms)
    assert output["double_relax"][1]["parameters"]["encut"] == 520
    assert output["double_relax"][1]["parameters"]["isym"] == 0
    assert output["double_relax"][1]["parameters"]["nsw"] > 0
    assert output["double_relax"][1]["parameters"]["isif"] == 3

    assert output["nsites"] == len(atoms)
    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["sigma"] == 0.01
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["laechg"] is True

    output = qmof_relax_job(atoms, run_prerelax=False)
    assert output["prerelax_lowacc"] is None

    output = qmof_relax_job(atoms, preset="BulkSet", nelmin=6)
    assert output["double_relax"][0]["parameters"]["encut"] == 520
    assert output["double_relax"][0]["parameters"]["nelmin"] == 6
    assert output["double_relax"][0]["parameters"]["sigma"] == 0.05

    assert output["double_relax"][1]["parameters"]["encut"] == 520
    assert output["double_relax"][1]["parameters"]["nelmin"] == 6
    assert output["double_relax"][1]["parameters"]["sigma"] == 0.05

    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["nelmin"] == 6
    assert output["parameters"]["sigma"] == 0.05

    output = qmof_relax_job(atoms, relax_cell=False)
    assert "volume-relax" not in output

    assert output["double_relax"][0]["parameters"]["isif"] == 2
    assert output["double_relax"][1]["parameters"]["isif"] == 2


def test_mp_metagga_prerelax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Al")
    output = mp_metagga_prerelax_job(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"] == {
        "algo": "all",
        "ediff": 1e-5,
        "ediffg": -0.05,
        "efermi": "midgap",  # added by copilot
        "enaug": 1360,
        "encut": 680,
        "gga": "ps",
        "ibrion": 2,
        "isif": 3,
        "ismear": 0,
        "ispin": 2,
        "kspacing": 0.22,
        "laechg": False,  # disabled by us
        "lasph": True,
        "lcharg": True,
        "lelf": False,
        "lmixtau": True,
        "lorbit": 11,
        "lreal": "auto",
        "lvtot": False,  # disabled by us
        "lwave": True,
        "magmom": [0.6],
        "nelm": 200,
        "nsw": 99,
        "prec": "accurate",
        "setups": {"Al": ""},
        "sigma": 0.05,
        "pp": "pbe",
    }

    output = mp_metagga_prerelax_job(atoms, bandgap=0)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["gga"] == "ps"
    assert output["parameters"]["ediffg"] == -0.05
    assert output["parameters"]["encut"] == 680
    assert output["parameters"]["kspacing"] == 0.22
    assert output["parameters"]["ismear"] == 0
    assert output["parameters"]["sigma"] == 0.05
    assert output["parameters"]["pp"] == "pbe"
    assert "metagga" not in output["parameters"]

    output = mp_metagga_prerelax_job(atoms, bandgap=100)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["gga"] == "ps"
    assert output["parameters"]["ediffg"] == -0.05
    assert output["parameters"]["encut"] == 680
    assert output["parameters"]["kspacing"] == 0.44
    assert output["parameters"]["ismear"] == 0
    assert output["parameters"]["sigma"] == 0.05
    assert output["parameters"]["pp"] == "pbe"
    assert "metagga" not in output["parameters"]


def test_mp_metagga_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Al")

    output = mp_metagga_relax_job(atoms)
    assert output["relax2"]["nsites"] == len(atoms)
    assert output["relax2"]["parameters"] == {
        "algo": "all",
        "ediff": 1e-5,
        "ediffg": -0.02,
        "efermi": "midgap",  # added by copilot
        "enaug": 1360,
        "encut": 680,
        "ibrion": 2,
        "isif": 3,
        "ismear": 0,
        "ispin": 2,
        "kspacing": 0.22,
        "laechg": False,  # disabled by us
        "lasph": True,
        "lcharg": True,
        "lelf": False,
        "lmixtau": True,
        "lorbit": 11,
        "lreal": "auto",
        "lvtot": False,  # disabled by us
        "lwave": True,
        "magmom": [0.6],
        "metagga": "r2scan",
        "nelm": 200,
        "nsw": 99,
        "prec": "accurate",
        "sigma": 0.05,
        "pp": "pbe",
        "setups": {"Al": ""},
    }

    output = mp_metagga_relax_job(atoms, bandgap=0)
    assert output["relax2"]["nsites"] == len(atoms)
    assert output["relax2"]["parameters"]["metagga"].lower() == "r2scan"
    assert output["relax2"]["parameters"]["ediffg"] == -0.02
    assert output["relax2"]["parameters"]["encut"] == 680
    assert output["relax2"]["parameters"]["kspacing"] == 0.22
    assert output["relax2"]["parameters"]["ismear"] == 0
    assert output["relax2"]["parameters"]["sigma"] == 0.05
    assert output["relax2"]["parameters"]["pp"] == "pbe"

    output = mp_metagga_relax_job(atoms, bandgap=100)
    assert output["relax2"]["nsites"] == len(atoms)
    assert output["relax2"]["parameters"]["metagga"].lower() == "r2scan"
    assert output["relax2"]["parameters"]["ediffg"] == -0.02
    assert output["relax2"]["parameters"]["encut"] == 680
    assert output["relax2"]["parameters"]["kspacing"] == 0.44
    assert output["relax2"]["parameters"]["ismear"] == 0
    assert output["relax2"]["parameters"]["sigma"] == 0.05
    assert output["relax2"]["parameters"]["pp"] == "pbe"


def test_mp_metagga_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Al")

    output = mp_metagga_static_job(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"] == {
        "algo": "fast",
        "ediff": 1e-05,
        "efermi": "midgap",  # added by copilot
        "enaug": 1360,
        "encut": 680,
        "ismear": -5,
        "ispin": 2,
        "kspacing": 0.22,
        "laechg": True,
        "lasph": True,
        "lcharg": True,
        "lelf": False,
        "lmixtau": True,
        "lorbit": 11,
        "lreal": False,
        "lvtot": True,
        "lwave": True,  # enabled by us
        "magmom": [0.6],
        "metagga": "r2scan",
        "nelm": 200,
        "nsw": 0,
        "prec": "accurate",
        "sigma": 0.05,
        "pp": "pbe",
        "setups": {"Al": ""},
    }


def test_mp_metagga_relax_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Al")

    output = mp_metagga_relax_flow(atoms)
    assert output["static"]["nsites"] == len(atoms)
    assert output["relax"]["relax2"]["parameters"]["metagga"].lower() == "r2scan"
    assert output["relax"]["relax2"]["parameters"]["ediffg"] == -0.02
    assert output["relax"]["relax2"]["parameters"]["encut"] == 680
    assert output["relax"]["relax2"]["parameters"]["ismear"] == 0
    assert output["relax"]["relax2"]["parameters"]["sigma"] == 0.05
    assert output["relax"]["relax2"]["parameters"]["kspacing"] == 0.22
    assert output["relax"]["relax2"]["parameters"]["pp"] == "pbe"
    assert output["prerelax"]["parameters"]["gga"] == "ps"
    assert output["prerelax"]["parameters"]["ismear"] == 0
    assert output["prerelax"]["parameters"]["pp"] == "pbe"

    atoms = bulk("C")
    output = mp_metagga_relax_flow(atoms)
    assert output["static"]["nsites"] == len(atoms)
    assert output["relax"]["relax2"]["parameters"]["metagga"].lower() == "r2scan"
    assert output["relax"]["relax2"]["parameters"]["ediffg"] == -0.02
    assert output["relax"]["relax2"]["parameters"]["encut"] == 680
    assert output["relax"]["relax2"]["parameters"]["ismear"] == 0
    assert output["relax"]["relax2"]["parameters"]["kspacing"] == pytest.approx(
        0.28329488761304206
    )
    assert output["relax"]["relax2"]["parameters"]["pp"] == "pbe"
    assert output["prerelax"]["parameters"]["ismear"] == 0
    assert output["prerelax"]["parameters"]["pp"] == "pbe"

    atoms = molecule("O2")
    atoms.center(vacuum=10)
    atoms.pbc = True
    output = mp_metagga_relax_flow(atoms)
    assert output["static"]["nsites"] == len(atoms)
    assert output["static"]["parameters"]["ismear"] == -5
    assert output["static"]["parameters"]["nsw"] == 0
    assert output["static"]["parameters"]["algo"] == "fast"
    assert output["relax"]["relax2"]["parameters"]["metagga"].lower() == "r2scan"
    assert output["relax"]["relax2"]["parameters"]["ediffg"] == -0.02
    assert output["relax"]["relax2"]["parameters"]["encut"] == 680
    assert output["relax"]["relax2"]["parameters"]["ismear"] == 0
    assert output["relax"]["relax2"]["parameters"]["kspacing"] == pytest.approx(
        0.28329488761304206
    )
    assert output["relax"]["relax2"]["parameters"]["pp"] == "pbe"
    assert output["prerelax"]["parameters"]["ismear"] == 0
    assert output["prerelax"]["parameters"]["pp"] == "pbe"


def test_mp_gga_relax_job():
    atoms = bulk("Ni") * (2, 1, 1)
    atoms[0].symbol = "O"
    output = mp_gga_relax_job(atoms)
    assert output["relax2"]["nsites"] == len(atoms)
    assert output["relax1"]["atoms"].get_chemical_symbols() == ["O", "Ni"]
    assert output["relax2"]["atoms"].get_chemical_symbols() == ["O", "Ni"]
    assert output["relax2"]["parameters"] == {
        "algo": "fast",
        "ediff": 0.0001,
        "efermi": "midgap",  # added by copilot
        "encut": 520,
        "gamma": True,
        "ibrion": 2,
        "isif": 3,
        "ismear": -5,
        "ispin": 2,
        "kpts": [5, 11, 11],
        "lasph": True,
        "ldau": True,
        "ldauj": [0, 0],
        "ldaul": [0, 2],
        "ldauprint": 1,  # added by us (sensible)
        "ldautype": 2,  # added by us (sensible)
        "ldauu": [0, 6.2],
        "lmaxmix": 4,
        "lorbit": 11,
        "lreal": "auto",
        "lwave": False,
        "magmom": [0.6, 0.6],
        "nelm": 100,
        "nsw": 99,
        "prec": "accurate",
        "sigma": 0.05,
        "pp": "pbe",
        "setups": {"O": "", "Ni": "_pv"},
    }


def test_mp_gga_static_job():
    atoms = bulk("Ni") * (2, 1, 1)
    atoms[0].symbol = "O"
    output = mp_gga_static_job(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"] == {
        "algo": "fast",
        "ediff": 0.0001,
        "efermi": "midgap",  # added by copilot
        "encut": 520,
        "gamma": True,
        "ismear": -5,
        "ispin": 2,
        "kpts": [6, 13, 13],
        "lasph": True,
        "lcharg": True,  # modified by us (sensible)
        "ldau": True,
        "ldauj": [0, 0],
        "ldaul": [0, 2],
        "ldauprint": 1,  # added by us (sensible)
        "ldautype": 2,  # added by us (sensible)
        "ldauu": [0, 6.2],
        "lmaxmix": 4,
        "lorbit": 11,
        "lreal": False,
        "lwave": True,  # modified by us (sensible)
        "magmom": [0.6, 0.6],
        "nelm": 100,
        "nsw": 0,
        "prec": "accurate",
        "sigma": 0.05,
        "pp": "pbe",
        "setups": {"Ni": "_pv", "O": ""},
    }


def test_mp_gga_relax_flow():
    atoms = bulk("Ni") * (2, 1, 1)
    atoms[0].symbol = "O"
    output = mp_gga_relax_flow(atoms)
    assert output["static"]["nsites"] == len(atoms)
    relax_params = {
        "algo": "fast",
        "ediff": 0.0001,
        "efermi": "midgap",  # added by copilot
        "encut": 520,
        "gamma": True,
        "ibrion": 2,
        "isif": 3,
        "ismear": -5,
        "ispin": 2,
        "kpts": [5, 11, 11],
        "lasph": True,
        "ldau": True,
        "ldauj": [0, 0],
        "ldaul": [0, 2],
        "ldauprint": 1,  # added by us (sensible)
        "ldautype": 2,  # added by us (sensible)
        "ldauu": [0, 6.2],
        "lmaxmix": 4,
        "lorbit": 11,
        "lreal": "auto",
        "lwave": False,
        "magmom": [0.6, 0.6],
        "nelm": 100,
        "nsw": 99,
        "prec": "accurate",
        "sigma": 0.05,
        "pp": "pbe",
        "setups": {"O": "", "Ni": "_pv"},
    }
    assert output["relax"]["relax1"]["parameters"] == relax_params
    assert output["relax"]["relax2"]["parameters"] == relax_params
    assert output["static"]["parameters"] == {
        "algo": "fast",
        "ediff": 0.0001,
        "efermi": "midgap",  # added by copilot
        "encut": 520,
        "gamma": True,
        "ismear": -5,
        "ispin": 2,
        "kpts": [6, 13, 13],
        "lasph": True,
        "lcharg": True,  # modified by us (sensible)
        "ldau": True,
        "ldauj": [0, 0],
        "ldaul": [0, 2],
        "ldauprint": 1,  # added by us (sensible)
        "ldautype": 2,  # added by us (sensible)
        "ldauu": [0, 6.2],
        "lmaxmix": 4,
        "lorbit": 11,
        "lreal": False,
        "lwave": True,  # modified by us (sensible)
        "magmom": [0.6, 0.6],
        "nelm": 100,
        "nsw": 0,
        "prec": "accurate",
        "sigma": 0.05,
        "pp": "pbe",
        "setups": {"Ni": "_pv", "O": ""},
    }


def test_mp_gga_relax_flow_custom():
    atoms = bulk("Ni") * (2, 1, 1)
    atoms[0].symbol = "O"
    output = mp_metagga_relax_flow(
        mp_gga_relax_flow(
            atoms["static"]["atoms"], job_params={"mp_gga_relax_job": {"nsw": 0}}
        ),
        job_params={"mp_metagga_relax_job": {"nsw": 0}},
    )
    assert output["relax"]["parameters"]["nsw"] == 0
