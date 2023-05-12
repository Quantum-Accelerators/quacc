import os
from shutil import rmtree

import pytest
from ase.build import bulk, molecule

from quacc.recipes.vasp.core import double_relax_job, relax_job, static_job
from quacc.recipes.vasp.qmof import qmof_relax_job
from quacc.recipes.vasp.slabs import (
    BulkToSlabsFlow,
    SlabToAdsFlow,
    slab_relax_job,
    slab_static_job,
)


def teardown_module():
    for f in os.listdir(os.getcwd()):
        if "quacc-tmp" in f or f == "tmp_dir":
            if os.path.islink(f):
                os.unlink(f)
            else:
                rmtree(f)


def test_static_job():
    atoms = bulk("Cu") * (2, 2, 2)

    output = static_job(atoms)
    assert output["nsites"] == len(atoms)
    assert "isym" not in output["parameters"]
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["lwave"] == True

    output = static_job(atoms, preset="BulkSet", swaps={"ncore": 2, "kpar": 4})
    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["ncore"] == 2
    assert output["parameters"]["kpar"] == 4

    output = static_job(
        atoms, preset="QMOFSet", swaps={"ismear": 0, "sigma": 0.01, "nedos": None}
    )
    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["ismear"] == 0
    assert output["parameters"]["sigma"] == 0.01

    output = static_job(atoms, swaps={"lwave": None})
    assert "lwave" not in output["parameters"]


def test_relax_job():
    atoms = bulk("Cu") * (2, 2, 2)

    output = relax_job(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isif"] == 3
    assert output["parameters"]["lwave"] == False

    output = relax_job(atoms, preset="BulkSet", swaps={"nelmin": 6})
    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["nelmin"] == 6

    output = relax_job(atoms, volume_relax=False)
    assert output["parameters"]["isif"] == 2


def test_doublerelax_job():
    atoms = bulk("Cu") * (2, 2, 2)

    output = double_relax_job(atoms)
    assert output["relax1"]["nsites"] == len(atoms)
    assert output["relax1"]["parameters"]["isym"] == 0
    assert output["relax1"]["parameters"]["nsw"] > 0
    assert output["relax1"]["parameters"]["isif"] == 3
    assert output["relax1"]["parameters"]["lwave"] == True
    assert output["relax2"]["nsites"] == len(atoms)
    assert output["relax2"]["parameters"]["isym"] == 0
    assert output["relax2"]["parameters"]["nsw"] > 0
    assert output["relax2"]["parameters"]["isif"] == 3
    assert output["relax2"]["parameters"]["lwave"] == True

    output = double_relax_job(atoms, preset="BulkSet", swaps2={"nelmin": 6})
    assert output["relax1"]["parameters"]["encut"] == 520
    assert "nelmin" not in output["relax1"]["parameters"]
    assert output["relax2"]["parameters"]["encut"] == 520
    assert output["relax2"]["parameters"]["nelmin"] == 6

    output = double_relax_job(atoms, volume_relax=False)
    assert output["relax1"]["parameters"]["isif"] == 2
    assert output["relax2"]["parameters"]["isif"] == 2

    output = double_relax_job(atoms, swaps1={"kpts": [1, 1, 1]})


def test_slab_static_job():
    atoms = bulk("Cu") * (2, 2, 2)

    output = slab_static_job(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["idipol"] == 3
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["lvhar"] == True

    output = slab_static_job(atoms, preset="SlabSet", swaps={"nelmin": 6})
    assert output["parameters"]["encut"] == 450
    assert output["parameters"]["nelmin"] == 6

    output = slab_static_job(atoms, preset="SlabSet", swaps={"encut": None})
    assert "encut" not in output["parameters"]


def test_slab_relax_job():
    atoms = bulk("Cu") * (2, 2, 2)

    output = slab_relax_job(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isif"] == 2
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["lwave"] == False

    output = slab_relax_job(atoms, preset="SlabSet", swaps={"nelmin": 6})
    assert output["parameters"]["encut"] == 450
    assert output["parameters"]["nelmin"] == 6


def test_slab_dynamic_jobs():
    atoms = bulk("Cu") * (2, 2, 2)

    ### --------- Test BulkToSlabsFlow --------- ###

    with pytest.raises(ValueError):
        BulkToSlabsFlow(relax_electron=None, static_electron=None).run(atoms)

    outputs = BulkToSlabsFlow(static_electron=None).run(atoms)
    assert outputs[0]["nsites"] == 96
    assert [output["parameters"]["isif"] == 2 for output in outputs]

    outputs = BulkToSlabsFlow().run(atoms)
    assert outputs[0]["nsites"] == 96
    assert [output["parameters"]["nsw"] == 0 for output in outputs]

    outputs = BulkToSlabsFlow(
        relax_kwargs={"preset": "SlabSet", "swaps": {"nelmin": 6}},
        static_electron=None,
    ).run(atoms)

    assert outputs[0]["nsites"] == 96
    assert [output["parameters"]["isif"] == 2 for output in outputs]
    assert [output["parameters"]["nelmin"] == 6 for output in outputs]
    assert [output["parameters"]["encut"] == 450 for output in outputs]

    outputs = BulkToSlabsFlow(
        relax_electron=None,
        static_kwargs={"preset": "SlabSet", "swaps": {"nelmin": 6}},
    ).run(atoms)

    assert outputs[0]["nsites"] == 96
    assert [output["parameters"]["nsw"] == 0 for output in outputs]
    assert [output["parameters"]["nelmin"] == 6 for output in outputs]
    assert [output["parameters"]["encut"] == 450 for output in outputs]

    outputs = BulkToSlabsFlow(
        static_kwargs={"preset": "SlabSet", "swaps": {"nelmin": 6}},
    ).run(atoms)

    assert outputs[0]["nsites"] == 96
    assert [output["parameters"]["nsw"] == 0 for output in outputs]
    assert [output["parameters"]["nelmin"] == 6 for output in outputs]
    assert [output["parameters"]["encut"] == 450 for output in outputs]

    ### --------- Test SlabToAdsorbatesJob --------- ###
    atoms = outputs[0]["atoms"]
    adsorbate = molecule("H2")

    outputs = SlabToAdsFlow(static_electron=None).run(atoms, adsorbate)

    assert outputs[0]["nsites"] == 98
    assert [output["parameters"]["isif"] == 2 for output in outputs]

    outputs = SlabToAdsFlow().run(atoms, adsorbate)
    assert outputs[0]["nsites"] == 98
    assert [output["parameters"]["nsw"] == 0 for output in outputs]

    outputs = SlabToAdsFlow(
        relax_kwargs={"preset": "SlabSet", "swaps": {"nelmin": 6}},
        static_electron=None,
    ).run(atoms, adsorbate)

    assert outputs[0]["nsites"] == 98
    assert [output["parameters"]["isif"] == 2 for output in outputs]
    assert [output["parameters"]["nelmin"] == 6 for output in outputs]
    assert [output["parameters"]["encut"] == 450 for output in outputs]

    outputs = SlabToAdsFlow(
        relax_electron=None,
        static_kwargs={"preset": "SlabSet", "swaps": {"nelmin": 6}},
    ).run(atoms, adsorbate)

    assert outputs[0]["nsites"] == 98
    assert [output["parameters"]["nsw"] == 0 for output in outputs]
    assert [output["parameters"]["nelmin"] == 6 for output in outputs]
    assert [output["parameters"]["encut"] == 450 for output in outputs]

    outputs = SlabToAdsFlow(
        static_kwargs={"preset": "SlabSet", "swaps": {"nelmin": 6}},
    ).run(atoms, adsorbate)

    assert outputs[0]["nsites"] == 98
    assert [output["parameters"]["nsw"] == 0 for output in outputs]
    assert [output["parameters"]["nelmin"] == 6 for output in outputs]
    assert [output["parameters"]["encut"] == 450 for output in outputs]

    adsorbate2 = molecule("CH3")
    adsorbate2.set_initial_magnetic_moments([1, 0, 0, 0])
    outputs = SlabToAdsFlow().run(atoms, adsorbate2)
    assert outputs[0]["nsites"] == 99
    assert [output["parameters"]["nsw"] == 0 for output in outputs]

    outputs = SlabToAdsFlow(
        relax_kwargs={"preset": "SlabSet", "swaps": {"nelmin": 6}},
        static_kwargs={"preset": "SlabSet", "swaps": {"nelmin": 6}},
    ).run(atoms, adsorbate)

    assert outputs[0]["nsites"] == 98
    assert [output["parameters"]["nsw"] == 0 for output in outputs]
    assert [output["parameters"]["nelmin"] == 6 for output in outputs]
    assert [output["parameters"]["encut"] == 450 for output in outputs]


def test_qmof():
    atoms = bulk("Cu")
    output = qmof_relax_job(atoms)
    assert output["prerelax-lowacc"]["nsites"] == len(atoms)
    assert output["prerelax-lowacc"]["parameters"]["sigma"] == 0.01
    assert output["prerelax-lowacc"]["parameters"]["isym"] == 0
    assert output["prerelax-lowacc"]["parameters"]["nsw"] == 0
    assert "isif" not in output["prerelax-lowacc"]["parameters"]
    assert "encut" not in output["prerelax-lowacc"]["parameters"]

    assert output["position-relax-lowacc"]["nsites"] == len(atoms)
    assert output["position-relax-lowacc"]["parameters"]["sigma"] == 0.01
    assert output["position-relax-lowacc"]["parameters"]["isym"] == 0
    assert output["position-relax-lowacc"]["parameters"]["nsw"] > 0
    assert output["position-relax-lowacc"]["parameters"]["isif"] == 2
    assert "encut" not in output["prerelax-lowacc"]["parameters"]

    assert output["volume-relax-lowacc"]["nsites"] == len(atoms)
    assert output["volume-relax-lowacc"]["parameters"]["encut"] == 520
    assert output["volume-relax-lowacc"]["parameters"]["sigma"] == 0.01
    assert output["volume-relax-lowacc"]["parameters"]["isym"] == 0
    assert output["volume-relax-lowacc"]["parameters"]["nsw"] > 0
    assert output["volume-relax-lowacc"]["parameters"]["isif"] == 3

    assert output["double-relax"]["relax1"]["nsites"] == len(atoms)
    assert output["double-relax"]["relax1"]["parameters"]["encut"] == 520
    assert output["double-relax"]["relax1"]["parameters"]["sigma"] == 0.01
    assert output["double-relax"]["relax1"]["parameters"]["isym"] == 0
    assert output["double-relax"]["relax1"]["parameters"]["nsw"] > 0
    assert output["double-relax"]["relax1"]["parameters"]["isif"] == 3

    assert output["double-relax"]["relax2"]["nsites"] == len(atoms)
    assert output["double-relax"]["relax2"]["parameters"]["encut"] == 520
    assert output["double-relax"]["relax2"]["parameters"]["isym"] == 0
    assert output["double-relax"]["relax2"]["parameters"]["nsw"] > 0
    assert output["double-relax"]["relax2"]["parameters"]["isif"] == 3

    assert output["static"]["nsites"] == len(atoms)
    assert output["static"]["parameters"]["encut"] == 520
    assert output["static"]["parameters"]["sigma"] == 0.01
    assert output["static"]["parameters"]["isym"] == 0
    assert output["static"]["parameters"]["nsw"] == 0
    assert output["static"]["parameters"]["laechg"] == True

    output = qmof_relax_job(atoms, prerelax=False)
    assert output["prerelax-lowacc"] is None

    output = qmof_relax_job(atoms, preset="BulkSet", swaps={"nelmin": 6})
    assert output["double-relax"]["relax1"]["parameters"]["encut"] == 520
    assert output["double-relax"]["relax1"]["parameters"]["nelmin"] == 6
    assert output["double-relax"]["relax1"]["parameters"]["sigma"] == 0.05

    assert output["double-relax"]["relax2"]["parameters"]["encut"] == 520
    assert output["double-relax"]["relax2"]["parameters"]["nelmin"] == 6
    assert output["double-relax"]["relax2"]["parameters"]["sigma"] == 0.05

    assert output["static"]["parameters"]["encut"] == 520
    assert output["static"]["parameters"]["nelmin"] == 6
    assert output["static"]["parameters"]["sigma"] == 0.05

    output = qmof_relax_job(atoms, volume_relax=False)
    assert "volume-relax" not in output

    assert output["double-relax"]["relax1"]["parameters"]["isif"] == 2
    assert output["double-relax"]["relax2"]["parameters"]["isif"] == 2

    atoms = bulk("Cu") * (8, 8, 8)
    output = qmof_relax_job(atoms)
