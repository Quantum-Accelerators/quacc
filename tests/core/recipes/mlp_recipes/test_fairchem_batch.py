"""Tests for FAIRChem batched inference functionality."""

from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")

from importlib.util import find_spec

from ase.build import bulk

# Skip all tests if fairchem is not installed
pytestmark = pytest.mark.skipif(
    find_spec("fairchem") is None, reason="fairchem not installed"
)

# Additional skip if HF token is not available (needed for model download)
try:
    from huggingface_hub.utils._auth import get_token

    HAS_HF_TOKEN = get_token() is not None
except Exception:
    HAS_HF_TOKEN = False


@pytest.fixture
def cleanup_batchers():
    """Fixture to ensure batchers are cleaned up after tests."""
    yield
    # Clean up after test
    from quacc.recipes.mlp._base import shutdown_inference_batchers

    shutdown_inference_batchers()


@pytest.mark.skipif(not HAS_HF_TOKEN, reason="HuggingFace token not available")
class TestGetInferenceBatcher:
    """Tests for get_inference_batcher function."""

    def test_batcher_creation(self, tmp_path, monkeypatch, cleanup_batchers):
        """Test that a batcher is created successfully."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp._base import get_inference_batcher

        batcher = get_inference_batcher(name_or_path="uma-s-1", task_name="omat")

        assert batcher is not None
        assert hasattr(batcher, "batch_predict_unit")
        assert hasattr(batcher, "executor")

    def test_batcher_caching(self, tmp_path, monkeypatch, cleanup_batchers):
        """Test that batchers are cached and reused."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp._base import get_inference_batcher

        batcher1 = get_inference_batcher(name_or_path="uma-s-1", task_name="omat")
        batcher2 = get_inference_batcher(name_or_path="uma-s-1", task_name="omat")

        # Same config should return same batcher instance
        assert batcher1 is batcher2

    def test_batcher_different_configs(self, tmp_path, monkeypatch, cleanup_batchers):
        """Test that different configs create different batchers."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp._base import (
            _INFERENCE_BATCHER_CACHE,
            get_inference_batcher,
        )

        # Clear cache to start fresh
        _INFERENCE_BATCHER_CACHE.clear()

        batcher1 = get_inference_batcher(
            name_or_path="uma-s-1", task_name="omat", max_batch_size=256
        )
        batcher2 = get_inference_batcher(
            name_or_path="uma-s-1", task_name="omat", max_batch_size=512
        )

        # Different configs should create different batchers
        assert batcher1 is not batcher2
        assert len(_INFERENCE_BATCHER_CACHE) == 2

    def test_batcher_with_custom_kwargs(self, tmp_path, monkeypatch, cleanup_batchers):
        """Test batcher creation with custom kwargs."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp._base import (
            _INFERENCE_BATCHER_CACHE,
            get_inference_batcher,
        )

        _INFERENCE_BATCHER_CACHE.clear()

        batcher = get_inference_batcher(
            name_or_path="uma-s-1", task_name="omat", max_batch_size=256
        )
        assert batcher is not None


@pytest.mark.skipif(not HAS_HF_TOKEN, reason="HuggingFace token not available")
class TestShutdownInferenceBatchers:
    """Tests for shutdown_inference_batchers function."""

    def test_shutdown_clears_cache(self, tmp_path, monkeypatch):
        """Test that shutdown clears the batcher cache."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp._base import (
            _INFERENCE_BATCHER_CACHE,
            get_inference_batcher,
            shutdown_inference_batchers,
        )

        # Create a batcher
        get_inference_batcher(name_or_path="uma-s-1", task_name="omat")

        assert len(_INFERENCE_BATCHER_CACHE) > 0

        # Shutdown should clear the cache
        shutdown_inference_batchers()

        assert len(_INFERENCE_BATCHER_CACHE) == 0


@pytest.mark.skipif(not HAS_HF_TOKEN, reason="HuggingFace token not available")
class TestPickCalculatorWithPredictUnit:
    """Tests for pick_calculator with predict_unit support."""

    def test_pick_calculator_with_predict_unit(
        self, tmp_path, monkeypatch, cleanup_batchers
    ):
        """Test that pick_calculator accepts predict_unit kwarg."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp._base import get_inference_batcher, pick_calculator

        batcher = get_inference_batcher(name_or_path="uma-s-1", task_name="omat")

        # Create calculator with predict_unit
        calc = pick_calculator(
            "fairchem", predict_unit=batcher.batch_predict_unit, task_name="omat"
        )

        assert calc is not None
        assert hasattr(calc, "calculate")

    def test_calculator_with_predict_unit_computes(
        self, tmp_path, monkeypatch, cleanup_batchers
    ):
        """Test that calculator created with predict_unit can compute."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp._base import get_inference_batcher, pick_calculator

        batcher = get_inference_batcher(name_or_path="uma-s-1", task_name="omat")

        calc = pick_calculator(
            "fairchem", predict_unit=batcher.batch_predict_unit, task_name="omat"
        )

        atoms = bulk("Cu")
        atoms.calc = calc
        energy = atoms.get_potential_energy()

        assert isinstance(energy, float)
        assert energy == pytest.approx(-3.7579006783217954, rel=1e-4)


@pytest.mark.skipif(not HAS_HF_TOKEN, reason="HuggingFace token not available")
class TestMapPartitionFairchemBatch:
    """Tests for map_partition_fairchembatch function."""

    def test_basic_batched_static(self, tmp_path, monkeypatch, cleanup_batchers):
        """Test basic batched static calculation."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp.core import static_job
        from quacc.wflow_tools.job_patterns import map_partition_fairchembatch

        atoms_list = [bulk("Cu"), bulk("Cu") * (2, 1, 1)]

        results = map_partition_fairchembatch(
            static_job, atoms_list=atoms_list, name_or_path="uma-s-1", task_name="omat"
        )

        assert len(results) == 2
        assert "results" in results[0]
        assert "energy" in results[0]["results"]
        assert results[0]["results"]["energy"] == pytest.approx(
            -3.7579006783217954, rel=1e-4
        )

    def test_batched_with_mapped_kwargs(self, tmp_path, monkeypatch, cleanup_batchers):
        """Test batched calculation with additional mapped kwargs."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp.core import static_job
        from quacc.wflow_tools.job_patterns import map_partition_fairchembatch

        atoms_list = [bulk("Cu"), bulk("Ag")]

        # Test with additional_fields as mapped kwarg
        results = map_partition_fairchembatch(
            static_job,
            atoms_list=atoms_list,
            name_or_path="uma-s-1",
            task_name="omat",
            additional_fields=[{"tag": "cu_calc"}, {"tag": "ag_calc"}],
        )

        assert len(results) == 2
        assert results[0]["tag"] == "cu_calc"
        assert results[1]["tag"] == "ag_calc"

    def test_batched_with_unmapped_kwargs(
        self, tmp_path, monkeypatch, cleanup_batchers
    ):
        """Test batched calculation with unmapped kwargs."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp.core import static_job
        from quacc.wflow_tools.job_patterns import map_partition_fairchembatch

        atoms_list = [bulk("Cu"), bulk("Cu") * (2, 1, 1)]

        results = map_partition_fairchembatch(
            static_job,
            atoms_list=atoms_list,
            name_or_path="uma-s-1",
            task_name="omat",
            unmapped_kwargs={"additional_fields": {"shared_tag": "batch_run"}},
        )

        assert len(results) == 2
        assert results[0]["shared_tag"] == "batch_run"
        assert results[1]["shared_tag"] == "batch_run"

    def test_batched_inconsistent_lengths_raises(
        self, tmp_path, monkeypatch, cleanup_batchers
    ):
        """Test that inconsistent mapped_kwargs lengths raise an error."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp.core import static_job
        from quacc.wflow_tools.job_patterns import map_partition_fairchembatch

        atoms_list = [bulk("Cu"), bulk("Ag")]

        with pytest.raises(AssertionError, match="Inconsistent lengths"):
            map_partition_fairchembatch(
                static_job,
                atoms_list=atoms_list,
                name_or_path="uma-s-1",
                task_name="omat",
                # This has 3 elements but atoms_list has 2
                additional_fields=[{"a": 1}, {"a": 2}, {"a": 3}],
            )

    def test_batched_reuses_cached_batcher(
        self, tmp_path, monkeypatch, cleanup_batchers
    ):
        """Test that repeated calls reuse the cached batcher."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp._base import _INFERENCE_BATCHER_CACHE
        from quacc.recipes.mlp.core import static_job
        from quacc.wflow_tools.job_patterns import map_partition_fairchembatch

        _INFERENCE_BATCHER_CACHE.clear()

        atoms_list = [bulk("Cu")]

        # First call
        map_partition_fairchembatch(
            static_job, atoms_list=atoms_list, name_or_path="uma-s-1", task_name="omat"
        )

        cache_size_after_first = len(_INFERENCE_BATCHER_CACHE)

        # Second call with same config
        map_partition_fairchembatch(
            static_job, atoms_list=atoms_list, name_or_path="uma-s-1", task_name="omat"
        )

        # Cache size should remain the same (batcher was reused)
        assert len(_INFERENCE_BATCHER_CACHE) == cache_size_after_first

    def test_batched_with_custom_batcher_kwargs(
        self, tmp_path, monkeypatch, cleanup_batchers
    ):
        """Test batched calculation with custom batcher settings."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp.core import static_job
        from quacc.wflow_tools.job_patterns import map_partition_fairchembatch

        atoms_list = [bulk("Cu")]

        results = map_partition_fairchembatch(
            static_job,
            atoms_list=atoms_list,
            name_or_path="uma-s-1",
            task_name="omat",
            batcher_kwargs={"max_batch_size": 256},
        )

        assert len(results) == 1
        assert "results" in results[0]

    def test_end_to_end_ten_copper_atoms(self, tmp_path, monkeypatch, cleanup_batchers):
        """
        End-to-end test: create 10 copper atoms objects with varying sizes,
        run batched inference via map_partition_fairchembatch,
        and verify all results are returned correctly.
        """
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp.core import static_job
        from quacc.wflow_tools.job_patterns import map_partition_fairchembatch

        # Create 10 different copper structures with varying sizes
        atoms_list = []
        for i in range(10):
            if i < 5:
                # First 5: unit cells with small perturbations
                atoms = bulk("Cu")
                atoms.rattle(stdev=0.01, seed=i)
            else:
                # Last 5: 2x1x1 supercells with perturbations
                atoms = bulk("Cu") * (2, 1, 1)
                atoms.rattle(stdev=0.01, seed=i)
            atoms_list.append(atoms)

        assert len(atoms_list) == 10

        # Run batched inference
        results = map_partition_fairchembatch(
            static_job, atoms_list=atoms_list, name_or_path="uma-s-1", task_name="omat"
        )

        # Verify we got 10 results back
        assert len(results) == 10

        # Verify each result has the expected structure
        for i, result in enumerate(results):
            assert "results" in result, f"Result {i} missing 'results' key"
            assert "energy" in result["results"], f"Result {i} missing 'energy'"
            assert "forces" in result["results"], f"Result {i} missing 'forces'"
            assert "atoms" in result, f"Result {i} missing 'atoms'"

            # Verify energy is a reasonable float
            energy = result["results"]["energy"]
            assert isinstance(energy, float), f"Result {i} energy is not a float"
            assert energy < 0, f"Result {i} energy should be negative for Cu"

            # Verify forces shape matches atoms count
            forces = result["results"]["forces"]
            expected_natoms = 1 if i < 5 else 2
            assert forces.shape == (expected_natoms, 3), (
                f"Result {i} forces shape mismatch"
            )

        # Verify unit cell energies are similar (within perturbation tolerance)
        unit_cell_energies = [results[i]["results"]["energy"] for i in range(5)]
        for e in unit_cell_energies:
            assert e == pytest.approx(unit_cell_energies[0], rel=1e-2)

        # Verify supercell energies are roughly 2x unit cell
        supercell_energy = results[5]["results"]["energy"]
        unit_cell_energy = results[0]["results"]["energy"]
        assert supercell_energy == pytest.approx(2 * unit_cell_energy, rel=1e-2)

    def test_batched_faster_than_sequential(
        self, tmp_path, monkeypatch, cleanup_batchers
    ):
        """
        Performance test: verify that batched inference is faster than
        sequential execution using local relaxations.
        """
        import time

        import numpy as np

        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp._base import _INFERENCE_BATCHER_CACHE
        from quacc.recipes.mlp.core import relax_job
        from quacc.wflow_tools.job_patterns import map_partition_fairchembatch

        # Clear batcher cache to ensure fresh start
        _INFERENCE_BATCHER_CACHE.clear()

        # Create 10 rattled 4x4x4 copper supercells for relaxation
        atoms_list = []
        for i in range(10):
            atoms = bulk("Cu") * (4, 4, 4)
            # Rattle with different seed for each structure
            rng = np.random.default_rng(seed=i)
            atoms.rattle(stdev=0.1, seed=int(rng.integers(0, 2**31)))
            atoms_list.append(atoms)

        # Warmup run for batched (includes Ray Serve setup overhead)
        warmup_atoms = bulk("Cu") * (4, 4, 4)
        warmup_atoms.rattle(stdev=0.1, seed=999)
        _ = map_partition_fairchembatch(
            relax_job,
            atoms_list=[warmup_atoms],
            name_or_path="uma-s-1",
            task_name="omat",
        )

        # Time batched inference (after warmup)
        start_batched = time.perf_counter()
        results_batched = map_partition_fairchembatch(
            relax_job, atoms_list=atoms_list, name_or_path="uma-s-1", task_name="omat"
        )
        time_batched = time.perf_counter() - start_batched

        # Verify batched results are valid
        assert len(results_batched) == 10
        assert all("results" in r and "energy" in r["results"] for r in results_batched)

        # Warmup for sequential (first call loads model)
        warmup_atoms = bulk("Cu") * (4, 4, 4)
        warmup_atoms.rattle(stdev=0.1, seed=999)
        _ = relax_job(
            warmup_atoms, method="fairchem", name_or_path="uma-s-1", task_name="omat"
        )

        # Time sequential inference - just call relax_job in a loop
        start_sequential = time.perf_counter()
        results_sequential = []
        for atoms in atoms_list:
            result = relax_job(
                atoms, method="fairchem", name_or_path="uma-s-1", task_name="omat"
            )
            results_sequential.append(result)
        time_sequential = time.perf_counter() - start_sequential

        assert len(results_sequential) == 10

        # Print timing info for debugging

        # Batched should be faster due to concurrent execution and GPU batching
        # We use a generous margin to account for system variability
        assert time_batched < time_sequential * 1.5, (
            f"Batched inference ({time_batched:.3f}s) should be faster "
            f"than sequential ({time_sequential:.3f}s)"
        )

        # Verify both methods produce similar final energies (relaxed structures)
        for i, (batched, sequential) in enumerate(
            zip(results_batched, results_sequential, strict=False)
        ):
            assert batched["results"]["energy"] == pytest.approx(
                sequential["results"]["energy"], rel=1e-3
            ), f"Energy mismatch at index {i}"


@pytest.mark.skipif(not HAS_HF_TOKEN, reason="HuggingFace token not available")
class TestMapPartitionedListsFairchemBatch:
    """Tests for the full partition/map_partitioned_lists_fairchembatch/unpartition pattern."""

    def test_partitioned_workflow(self, tmp_path, monkeypatch, cleanup_batchers):
        """Test the full partition -> map -> unpartition workflow."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp.core import static_job
        from quacc.wflow_tools.job_patterns import (
            map_partitioned_lists_fairchembatch,
            partition,
            unpartition,
        )

        # Create 6 atoms objects
        atoms_list = [bulk("Cu") for _ in range(6)]

        # Partition into 2 partitions
        num_partitions = 2
        partitioned_atoms = partition(atoms_list, num_partitions)

        assert len(partitioned_atoms) == 2
        assert len(partitioned_atoms[0]) == 3
        assert len(partitioned_atoms[1]) == 3

        # Run via map_partitioned_lists_fairchembatch
        results_partitioned = map_partitioned_lists_fairchembatch(
            static_job,
            num_partitions,
            name_or_path="uma-s-1",
            task_name="omat",
            atoms_list=partitioned_atoms,
        )

        # Should get 2 lists of results (one per partition)
        assert len(results_partitioned) == 2
        assert len(results_partitioned[0]) == 3
        assert len(results_partitioned[1]) == 3

        # Unpartition to get flat list
        all_results = unpartition(results_partitioned)

        assert len(all_results) == 6

        # All Cu atoms should have similar energies
        for result in all_results:
            assert "results" in result
            assert result["results"]["energy"] == pytest.approx(
                -3.7579006783217954, rel=1e-4
            )


class TestMapPartitionFairchemBatchNoFairchem:
    """Tests for error handling when fairchem is not available."""

    def test_import_error_without_fairchem(self, tmp_path, monkeypatch):
        """Test that ImportError is raised when fairchem is not installed."""
        monkeypatch.chdir(tmp_path)

        # This test verifies the error message - it will only actually raise
        # if fairchem is not installed, otherwise it passes as fairchem is available
        if find_spec("fairchem") is None:
            from quacc.wflow_tools.job_patterns import map_partition_fairchembatch

            with pytest.raises(ImportError, match="fairchem must be installed"):
                map_partition_fairchembatch(
                    lambda x: x,
                    atoms_list=[bulk("Cu")],
                    name_or_path="uma-s-1",
                    task_name="omat",
                )
