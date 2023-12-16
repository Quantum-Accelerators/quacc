from functools import partial

import pytest
from ase.build import bulk

from quacc import SETTINGS

redun = pytest.importorskip("redun")
pytestmark = pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE != "redun",
    reason="This test requires the Redun workflow engine",
)


@pytest.fixture()
def scheduler():
    return redun.Scheduler()


from quacc.recipes.emt.core import relax_job  # skipcq: PYL-C0412
from quacc.recipes.emt.slabs import bulk_to_slabs_flow  # skipcq: PYL-C0412


def test_redun_functools(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    result = scheduler.run(
        bulk_to_slabs_flow(
            atoms, slab_relax_job=partial(relax_job, opt_params={"fmax": 0.1})
        )
    )
    assert len(result) == 4
    assert "atoms" in result[-1]
    assert result[-1]["fmax"] == 0.1
