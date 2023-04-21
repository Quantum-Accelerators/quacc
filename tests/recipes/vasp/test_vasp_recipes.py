import os
from shutil import rmtree

import pytest
from ase.build import bulk, molecule

from quacc.recipes.vasp.core import DoubleRelaxJob, RelaxJob, StaticJob
from quacc.recipes.vasp.qmof import QMOFRelaxJob
from quacc.recipes.vasp.slabs import (  # BulkToAdsorbatesFlow,; BulkToSlabsJob,; SlabToAdsorbatesJob,
    SlabRelaxJob,
    SlabStaticJob,
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

    output = StaticJob(atoms)
    assert output["nsites"] == len(atoms)
    assert "isym" not in output["parameters"]
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["lwave"] == True

    output = StaticJob(atoms, preset="BulkSet", swaps={"ncore": 2, "kpar": 4})
    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["ncore"] == 2
    assert output["parameters"]["kpar"] == 4

    output = StaticJob(
        atoms, preset="QMOFSet", swaps={"ismear": 0, "sigma": 0.01, "nedos": None}
    )
    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["ismear"] == 0
    assert output["parameters"]["sigma"] == 0.01

    output = StaticJob(atoms, swaps={"lwave": None})
    assert "lwave" not in output["parameters"]


def test_relax_job():
    atoms = bulk("Cu") * (2, 2, 2)

    output = RelaxJob(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isif"] == 3
    assert output["parameters"]["lwave"] == False

    output = RelaxJob(atoms, preset="BulkSet", swaps={"nelmin": 6})
    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["nelmin"] == 6

    output = RelaxJob(atoms, volume_relax=False)
    assert output["parameters"]["isif"] == 2


def test_doublerelax_job():
    atoms = bulk("Cu") * (2, 2, 2)

    output = DoubleRelaxJob(atoms)
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

    output = DoubleRelaxJob(atoms, preset="BulkSet", swaps2={"nelmin": 6})
    assert output["relax1"]["parameters"]["encut"] == 520
    assert "nelmin" not in output["relax1"]["parameters"]
    assert output["relax2"]["parameters"]["encut"] == 520
    assert output["relax2"]["parameters"]["nelmin"] == 6

    output = DoubleRelaxJob(atoms, volume_relax=False)
    assert output["relax1"]["parameters"]["isif"] == 2
    assert output["relax2"]["parameters"]["isif"] == 2

    output = DoubleRelaxJob(atoms, swaps1={"kpts": [1, 1, 1]})


def test_slab_static_job():
    atoms = bulk("Cu") * (2, 2, 2)

    output = SlabStaticJob(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["idipol"] == 3
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["lvhar"] == True

    output = SlabStaticJob(atoms, preset="SlabSet", swaps={"nelmin": 6})
    assert output["parameters"]["encut"] == 450
    assert output["parameters"]["nelmin"] == 6

    output = SlabStaticJob(atoms, preset="SlabSet", swaps={"encut": None})
    assert "encut" not in output["parameters"]


def test_slab_relax_job():
    atoms = bulk("Cu") * (2, 2, 2)

    output = SlabRelaxJob(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isif"] == 2
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["lwave"] == False

    output = SlabRelaxJob(atoms, preset="SlabSet", swaps={"nelmin": 6})
    assert output["parameters"]["encut"] == 450
    assert output["parameters"]["nelmin"] == 6


# def test_slab_dynamic_jobs():
#     atoms = bulk("Cu") * (2, 2, 2)

#     ### --------- Test BulkToSlabsJob --------- ###
#     flow = BulkToSlabsJob(atoms).make(atoms)
#     responses = run_locally(flow, ensure_success=True)

#     assert len(responses) == 9
#     uuids = list(responses.keys())

#     # First job is a dummy job to make slabs and should have no output
#     output0 = responses[uuids[0]][1].output
#     assert "generated_slabs" in output0
#     assert len(output0["generated_slabs"][0]) == 96

#     output1 = responses[uuids[1]][1].output
#     assert output1["nsites"] > len(atoms)
#     assert output1["parameters"]["isif"] == 2

#     output2 = responses[uuids[-1]][1].output
#     assert output2["nsites"] == output1["nsites"]
#     assert output2["parameters"]["nsw"] == 0

#     # Now try with kwargs
#     flow = BulkToSlabsJob(
#         slab_relax_job=SlabRelaxJob(preset="SlabSet", swaps={"nelmin": 6}),
#         slab_static_job=SlabStaticJob(preset="SlabSet", swaps={"nelmin": 6}),
#     ).make(atoms)
#     responses = run_locally(flow, ensure_success=True)

#     assert len(responses) == 9
#     uuids = list(responses.keys())

#     output0 = responses[uuids[0]][1].output
#     assert "generated_slabs" in output0
#     assert len(output0["generated_slabs"][0]) == 96

#     output1 = responses[uuids[1]][1].output
#     assert output1["parameters"]["isif"] == 2
#     assert output1["parameters"]["nelmin"] == 6
#     assert output1["parameters"]["encut"] == 450

#     output2 = responses[uuids[-1]][1].output
#     assert output2["parameters"]["nsw"] == 0
#     assert output2["parameters"]["nelmin"] == 6
#     assert output2["parameters"]["encut"] == 450

#     ### --------- Test SlabToAdsorbatesJob --------- ###
#     atoms = output2["atoms"]
#     adsorbate = molecule("H2")

#     flow = SlabToAdsorbatesJob().make(atoms, adsorbate)
#     responses = run_locally(flow, ensure_success=True)

#     assert len(responses) == 9
#     uuids = list(responses.keys())

#     # First job is a dummy job to make slabs and should have no output
#     output0 = responses[uuids[0]][1].output
#     assert "input_slabs" in output0
#     assert "generated_slab_ads" in output0
#     assert "H2" in output0["generated_slab_ads"]
#     assert len(output0["generated_slab_ads"]["H2"][0]) == 98

#     output1 = responses[uuids[1]][1].output
#     assert output1["nsites"] == len(output2["atoms"]) + 2
#     assert output1["parameters"]["isif"] == 2

#     output2 = responses[uuids[-1]][1].output
#     assert output2["nsites"] == output1["nsites"]
#     assert output2["parameters"]["nsw"] == 0

#     # Now try with kwargs
#     flow = SlabToAdsorbatesJob(
#         slab_ads_relax_job=SlabRelaxJob(preset="SlabSet", swaps={"nelmin": 6}),
#         slab_ads_static_job=SlabStaticJob(preset="SlabSet", swaps={"nelmin": 6}),
#     ).make(atoms, adsorbate)
#     responses = run_locally(flow, ensure_success=True)

#     assert len(responses) == 9
#     uuids = list(responses.keys())

#     output0 = responses[uuids[0]][1].output
#     assert "generated_slab_ads" in output0
#     assert "H2" in output0["generated_slab_ads"]
#     assert len(output0["generated_slab_ads"]["H2"][0]) == 98

#     output1 = responses[uuids[1]][1].output
#     assert output1["parameters"]["isif"] == 2
#     assert output1["parameters"]["nelmin"] == 6
#     assert output1["parameters"]["encut"] == 450

#     output2 = responses[uuids[-1]][1].output
#     assert output2["parameters"]["nsw"] == 0
#     assert output2["parameters"]["nelmin"] == 6
#     assert output2["parameters"]["encut"] == 450

#     # Now try with different adsorbate
#     adsorbate2 = molecule("CH3")
#     adsorbate2.set_initial_magnetic_moments([1, 0, 0, 0])
#     flow = SlabToAdsorbatesJob().make(atoms, adsorbate2)
#     responses = run_locally(flow, ensure_success=True)
#     assert len(responses) == 9

#     adsorbate2 = molecule("CH3")
#     flow = SlabToAdsorbatesJob().make(atoms, adsorbate2)
#     responses = run_locally(flow, ensure_success=True)
#     assert len(responses) == 9


# def test_slab_flows():
#     # TODO: This could use some more detailed tests

#     atoms = bulk("Cu")
#     adsorbate = molecule("H2O")
#     flow = BulkToAdsorbatesFlow().make(atoms, adsorbate, n_stable_slabs=1)
#     responses = run_locally(flow, ensure_success=True)
#     uuids = list(responses.keys())
#     assert len(responses) == 21

#     output0 = responses[uuids[0]][1].output
#     assert output0["parameters"]["ediffg"] == -0.02

#     output1 = responses[uuids[1]][1].output
#     assert output1["parameters"]["nsw"] == 0

#     output_final = responses[uuids[-1]][1].output
#     assert output_final["parameters"]["nsw"] == 0

#     flow = BulkToAdsorbatesFlow().make(
#         atoms, [adsorbate, molecule("H2")], n_stable_slabs=1
#     )
#     responses = run_locally(flow, ensure_success=True)
#     assert len(responses) == 29

#     flow = BulkToAdsorbatesFlow().make(atoms, adsorbate)
#     responses = run_locally(flow, ensure_success=True)
#     assert len(responses) == 48

#     flow = BulkToAdsorbatesFlow(bulk_relax_job=None, bulk_static_job=None).make(
#         atoms, adsorbate, n_stable_slabs=None
#     )
#     responses = run_locally(flow, ensure_success=True)
#     assert len(responses) == 46

#     with pytest.raises(ValueError):
#         flow = BulkToAdsorbatesFlow(bulk_relax_job=None, bulk_static_job=None).make(
#             atoms, adsorbate, n_stable_slabs=1
#         )


def test_qmof():
    atoms = bulk("Cu")
    output = QMOFRelaxJob(atoms)
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

    output = QMOFRelaxJob(atoms, prerelax=False)
    assert output["prerelax-lowacc"] is None

    output = QMOFRelaxJob(atoms, preset="BulkSet", swaps={"nelmin": 6})
    assert output["double-relax"]["relax1"]["parameters"]["encut"] == 520
    assert output["double-relax"]["relax1"]["parameters"]["nelmin"] == 6
    assert output["double-relax"]["relax1"]["parameters"]["sigma"] == 0.05

    assert output["double-relax"]["relax2"]["parameters"]["encut"] == 520
    assert output["double-relax"]["relax2"]["parameters"]["nelmin"] == 6
    assert output["double-relax"]["relax2"]["parameters"]["sigma"] == 0.05

    assert output["static"]["parameters"]["encut"] == 520
    assert output["static"]["parameters"]["nelmin"] == 6
    assert output["static"]["parameters"]["sigma"] == 0.05

    output = QMOFRelaxJob(atoms, volume_relax=False)
    assert "volume-relax" not in output

    assert output["double-relax"]["relax1"]["parameters"]["isif"] == 2
    assert output["double-relax"]["relax2"]["parameters"]["isif"] == 2

    atoms = bulk("Cu") * (8, 8, 8)
    output = QMOFRelaxJob(atoms)
