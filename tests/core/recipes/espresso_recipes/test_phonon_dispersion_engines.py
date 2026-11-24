from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.parametrize("engine", ["dask", "jobflow"])
def test_deferred_dispersion_flow(tmp_path, engine):
    """Resolve upstream atoms, path flags, file dependencies, and the final schema."""
    pytest.importorskip(engine)
    config = tmp_path / "quacc.yaml"
    config.write_text(f"WORKFLOW_ENGINE: {engine}\n")
    env = os.environ | {
        "QUACC_CONFIG_FILE": str(config),
        "QUACC_RESULTS_DIR": str(tmp_path),
        "QUACC_SCRATCH_DIR": str(tmp_path),
        "PYTHONPATH": str(Path(__file__).resolve().parents[4] / "src"),
    }
    script = r"""
import numpy as np
from ase import Atoms
from quacc import flow, job
from quacc.recipes.espresso import phonons

calls = []
def run_and_summarize(*, template, calc_defaults, calc_swaps, copy_files, **kwargs):
    binary = template.binary
    calls.append(binary)
    if binary == "q2r":
        assert copy_files == [{"source": "ph-results", "filenames": "matdyn*"}]
    if binary == "matdyn":
        assert copy_files == [{"source": "q2r-results", "filenames": "q2r.fc*"}]
        assert not calc_swaps["input_data"]["input"].get("q_in_band_form", False)
        assert np.asarray(calc_swaps["qpts"]).shape == (3, 3)
    results = {}
    if binary == "matdyn":
        results["matdyn_results"] = {"frequencies": [[-1, 2, 3], [4, 5, 6], [7, 8, 9]], "frequency_units": "cm^-1"}
    return {"dir_name": binary + "-results", "parameters": calc_swaps, "results": results}

phonons.run_and_summarize = run_and_summarize

@job
def scf():
    return {"atoms": Atoms("H", cell=[2, 3, 4], pbc=True), "dir_name": "scf-results"}

@flow
def workflow():
    previous = scf()
    return phonons.phonon_dispersion_flow(
        previous["atoms"], prev_outdir=previous["dir_name"],
        qpts=[[0, 0, 0, 2], [0.5, 0, 0, 1]], qpoint_labels=["G", "X"]
    )

from quacc import get_settings
if get_settings().WORKFLOW_ENGINE == "dask":
    result = workflow().compute(scheduler="synchronous")
else:
    from jobflow import JobStore, run_locally
    from maggma.stores import MemoryStore
    graph = workflow()
    responses = run_locally(graph, store=JobStore(MemoryStore()), ensure_success=True)
    outputs = [response.output for jobs in responses.values() for response in jobs.values()]
    result = next(output for output in outputs if isinstance(output, dict) and "dispersion" in output)
assert calls == ["ph", "q2r", "matdyn"], calls
assert np.asarray(result["dispersion"]["frequencies"]).shape == (3, 3)
assert result["dispersion"]["frequencies"][0][0] == -1
assert result["dispersion"]["tick_labels"] == ["G", "X"]
assert result["dispersion"]["metadata"]["qgrid"] == [4, 4, 4]
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
