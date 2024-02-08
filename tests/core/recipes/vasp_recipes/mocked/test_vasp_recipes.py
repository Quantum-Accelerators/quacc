import pytest
from ase.build import bulk, molecule

from quacc import SETTINGS
from quacc.recipes.vasp.core import double_relax_job, relax_job, static_job
from quacc.recipes.vasp.mp import (
    mp_metagga_prerelax_job,
    mp_metagga_relax_flow,
    mp_metagga_relax_job,
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
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isif"] == 3
    assert output["parameters"]["lwave"] is False
    assert output["relax1"]["parameters"]["encut"] == 520
    assert output["parameters"]["encut"] == 520

    output = double_relax_job(atoms, relax2_kwargs={"nelmin": 6})
    assert output["relax1"]["nsites"] == len(atoms)
    assert output["relax1"]["parameters"]["isym"] == 0
    assert output["relax1"]["parameters"]["nsw"] > 0
    assert output["relax1"]["parameters"]["isif"] == 3
    assert output["relax1"]["parameters"]["lwave"] is False
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isif"] == 3
    assert output["parameters"]["lwave"] is False
    assert output["relax1"]["parameters"]["encut"] == 520
    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["nelmin"] == 6

    output = double_relax_job(atoms, relax_cell=False)
    assert output["relax1"]["nsites"] == len(atoms)
    assert output["relax1"]["parameters"]["isym"] == 0
    assert output["relax1"]["parameters"]["nsw"] > 0
    assert output["relax1"]["parameters"]["lwave"] is False
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["lwave"] is False
    assert output["relax1"]["parameters"]["encut"] == 520
    assert output["parameters"]["encut"] == 520
    assert output["relax1"]["parameters"]["isif"] == 2
    assert output["parameters"]["isif"] == 2

    double_relax_job(atoms, relax1_kwargs={"kpts": [1, 1, 1]})


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
    assert output["parameters"]["gga"] == "PS"
    assert output["parameters"]["ediffg"] == -0.05
    assert output["parameters"]["encut"] == 680
    assert output["parameters"]["kspacing"] == 0.22
    assert output["parameters"]["ismear"] == 0
    assert output["parameters"]["sigma"] == 0.05
    assert output["parameters"]["pp"] == "PBE"
    assert "metagga" not in output["parameters"]

    output = mp_metagga_prerelax_job(atoms, bandgap=0)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["gga"] == "PS"
    assert output["parameters"]["ediffg"] == -0.05
    assert output["parameters"]["encut"] == 680
    assert output["parameters"]["kspacing"] == 0.22
    assert output["parameters"]["ismear"] == 2
    assert output["parameters"]["sigma"] == 0.2
    assert output["parameters"]["pp"] == "PBE"
    assert "metagga" not in output["parameters"]

    output = mp_metagga_prerelax_job(atoms, bandgap=100)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["gga"] == "PS"
    assert output["parameters"]["ediffg"] == -0.05
    assert output["parameters"]["encut"] == 680
    assert output["parameters"]["kspacing"] == 0.44
    assert output["parameters"]["ismear"] == -5
    assert output["parameters"]["sigma"] == 0.05
    assert output["parameters"]["pp"] == "PBE"
    assert "metagga" not in output["parameters"]


def test_mp_metagga_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Al")

    output = mp_metagga_relax_job(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["metagga"].lower() == "r2scan"
    assert output["parameters"]["ediffg"] == -0.02
    assert output["parameters"]["encut"] == 680
    assert output["parameters"]["kspacing"] == 0.22
    assert output["parameters"]["ismear"] == 0
    assert output["parameters"]["sigma"] == 0.05
    assert output["parameters"]["pp"] == "PBE"

    output = mp_metagga_relax_job(atoms, bandgap=0)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["metagga"].lower() == "r2scan"
    assert output["parameters"]["ediffg"] == -0.02
    assert output["parameters"]["encut"] == 680
    assert output["parameters"]["kspacing"] == 0.22
    assert output["parameters"]["ismear"] == 2
    assert output["parameters"]["sigma"] == 0.2
    assert output["parameters"]["pp"] == "PBE"

    output = mp_metagga_relax_job(atoms, bandgap=100)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["metagga"].lower() == "r2scan"
    assert output["parameters"]["ediffg"] == -0.02
    assert output["parameters"]["encut"] == 680
    assert output["parameters"]["kspacing"] == 0.44
    assert output["parameters"]["ismear"] == -5
    assert output["parameters"]["sigma"] == 0.05
    assert output["parameters"]["pp"] == "PBE"


def test_mp_metagga_relax_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Al")

    output = mp_metagga_relax_flow(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["metagga"].lower() == "r2scan"
    assert output["parameters"]["ediffg"] == -0.02
    assert output["parameters"]["encut"] == 680
    assert output["parameters"]["ismear"] == 2
    assert output["parameters"]["sigma"] == 0.2
    assert output["parameters"]["kspacing"] == 0.22
    assert output["parameters"]["pp"] == "PBE"
    assert output["prerelax"]["parameters"]["gga"] == "PS"
    assert output["prerelax"]["parameters"]["ismear"] == 0
    assert output["prerelax"]["parameters"]["pp"] == "PBE"

    atoms = bulk("C")
    output = mp_metagga_relax_flow(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["metagga"].lower() == "r2scan"
    assert output["parameters"]["ediffg"] == -0.02
    assert output["parameters"]["encut"] == 680
    assert output["parameters"]["ismear"] == -5
    assert output["parameters"]["kspacing"] == pytest.approx(0.28329488761304206)
    assert output["parameters"]["pp"] == "PBE"
    assert output["prerelax"]["parameters"]["ismear"] == 0
    assert output["prerelax"]["parameters"]["pp"] == "PBE"

    atoms = molecule("O2")
    atoms.center(vacuum=10)
    atoms.pbc = True
    output = mp_metagga_relax_flow(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["metagga"].lower() == "r2scan"
    assert output["parameters"]["ediffg"] == -0.02
    assert output["parameters"]["encut"] == 680
    assert output["parameters"]["ismear"] == -5
    assert output["parameters"]["kspacing"] == pytest.approx(0.28329488761304206)
    assert output["parameters"]["pp"] == "PBE"
    assert output["prerelax"]["parameters"]["ismear"] == 0
    assert output["prerelax"]["parameters"]["pp"] == "PBE"
