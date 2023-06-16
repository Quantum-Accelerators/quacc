import os
from shutil import rmtree

import numpy as np
import pytest
from ase.build import bulk

from quacc.recipes.emt.core import relax_job, static_job
from quacc.recipes.emt.jobflow.slabs import BulkToSlabsFlow as JFBulkToSlabsFlow

try:
    from chgnet.model.dynamics import CHGNetCalculator as CHGNet
except ImportError:
    CHGNet = None


def teardown_module():
    for f in os.listdir(os.getcwd()):
        if (
            f.endswith(".log")
            or f.endswith(".pckl")
            or f.endswith(".traj")
            or f.endswith(".out")
            or ".gz" in f
        ):
            os.remove(f)
        if "quacc-tmp" in f or f == "tmp_dir":
            rmtree(f)


@pytest.mark.skipif(CHGNet is None, reason="CHGNet must be installed.")
def test_static_Job():
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]

    output = static_job(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is False
    assert output["results"]["energy"] == pytest.approx(0.07001766638245854)


@pytest.mark.skipif(CHGNet is None, reason="CHGNet must be installed.")
def test_relax_Job():
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]

    output = relax_job(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is False
    assert output["results"]["energy"] == pytest.approx(-0.04543069081693929)
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01
