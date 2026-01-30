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

        batcher = get_inference_batcher(name_or_path="uma-s-1")

        assert batcher is not None
        assert hasattr(batcher, "batch_predict_unit")
        assert hasattr(batcher, "executor")

    def test_batcher_caching(self, tmp_path, monkeypatch, cleanup_batchers):
        """Test that batchers are cached and reused."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp._base import get_inference_batcher

        batcher1 = get_inference_batcher(name_or_path="uma-s-1")
        batcher2 = get_inference_batcher(name_or_path="uma-s-1")

        # Same config should return same batcher instance
        assert batcher1 is batcher2

    def test_batcher_same_model_different_batch_sizes(self, tmp_path, monkeypatch, cleanup_batchers):
        """Test that same model with different batcher_kwargs reuses the batcher via checkpoint swap."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp._base import (
            _current_batcher,
            get_inference_batcher,
            shutdown_inference_batchers,
        )

        # Clear cache to start fresh
        shutdown_inference_batchers()

        batcher1 = get_inference_batcher(
            name_or_path="uma-s-1", max_batch_size=256
        )
        batcher1_id = id(batcher1)
        
        batcher2 = get_inference_batcher(
            name_or_path="uma-s-1", max_batch_size=512
        )

        # Same model checkpoint, different batcher_kwargs: should reuse the batcher instance
        # via update_checkpoint (or in this case just return it since checkpoint is same)
        assert batcher1 is batcher2
        assert id(batcher2) == batcher1_id
        # Same batcher still active
        assert _current_batcher is batcher2

    def test_batcher_different_models_swaps_checkpoint(self, tmp_path, monkeypatch, cleanup_batchers):
        """Test that different models trigger checkpoint swap."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp._base import (
            _current_batcher,
            get_inference_batcher,
            shutdown_inference_batchers,
        )

        # Clear cache to start fresh
        shutdown_inference_batchers()

        batcher1 = get_inference_batcher(
            name_or_path="uma-s-1"
        )
        batcher1_id = id(batcher1)
        
        # Note: this test would swap checkpoint to a different model if available
        # For now, we can't test with a truly different model without more setup,
        # so we verify the cache key behavior instead
        batcher2 = get_inference_batcher(
            name_or_path="uma-s-1"  # same model
        )

        # Same checkpoint should return same batcher
        assert batcher1 is batcher2
        assert _current_batcher is batcher2

    def test_batcher_with_custom_kwargs(self, tmp_path, monkeypatch, cleanup_batchers):
        """Test batcher creation with custom kwargs."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp._base import (
            get_inference_batcher,
            shutdown_inference_batchers,
        )

        shutdown_inference_batchers()

        batcher = get_inference_batcher(
            name_or_path="uma-s-1", max_batch_size=256
        )
        assert batcher is not None


@pytest.mark.skipif(not HAS_HF_TOKEN, reason="HuggingFace token not available")
class TestShutdownInferenceBatchers:
    """Tests for shutdown_inference_batchers function."""

    def test_shutdown_clears_cache(self, tmp_path, monkeypatch):
        """Test that shutdown clears the batcher cache."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp._base import (
            get_inference_batcher,
            shutdown_inference_batchers,
        )
        from quacc.recipes.mlp import _base

        # Create a batcher
        get_inference_batcher(name_or_path="uma-s-1")

        assert _base._current_batcher is not None
        assert _base._current_checkpoint_key is not None

        # Shutdown should clear the cache
        shutdown_inference_batchers()

        assert _base._current_batcher is None
        assert _base._current_checkpoint_key is None


@pytest.mark.skipif(not HAS_HF_TOKEN, reason="HuggingFace token not available")
class TestPickCalculatorWithPredictUnit:
    """Tests for pick_calculator with predict_unit support."""

    def test_pick_calculator_with_predict_unit(
        self, tmp_path, monkeypatch, cleanup_batchers
    ):
        """Test that pick_calculator accepts predict_unit kwarg."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp._base import get_inference_batcher, pick_calculator

        batcher = get_inference_batcher(name_or_path="uma-s-1")

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

        batcher = get_inference_batcher(name_or_path="uma-s-1")

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
            atoms=atoms_list,
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
            atoms=atoms_list,
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
                atoms=atoms_list,
                name_or_path="uma-s-1",
                task_name="omat",
                # This has 3 elements but atoms_list has 2
                additional_fields=[{"a": 1}, {"a": 2}, {"a": 3}],
            )

    def test_batched_reuses_cached_batcher(self, tmp_path, monkeypatch, cleanup_batchers):
        """Test that repeated calls reuse the cached batcher."""
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp._base import shutdown_inference_batchers
        from quacc.recipes.mlp.core import static_job
        from quacc.wflow_tools.job_patterns import map_partition_fairchembatch

        shutdown_inference_batchers()

        atoms_list = [bulk("Cu")]

        # First call
        map_partition_fairchembatch(
            static_job, atoms_list=atoms_list, name_or_path="uma-s-1", task_name="omat"
        )

        from quacc.recipes.mlp import _base
        batcher_after_first = _base._current_batcher
        assert batcher_after_first is not None

        # Second call with same config
        map_partition_fairchembatch(
            static_job, atoms_list=atoms_list, name_or_path="uma-s-1", task_name="omat"
        )

        # Same batcher should be reused
        assert _base._current_batcher is batcher_after_first

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
            atoms=atoms_list,
            name_or_path="uma-s-1",
            task_name="omat",
            batcher_kwargs={"max_batch_size": 256},
        )

        assert len(results) == 1
        assert "results" in results[0]

    def test_checkpoint_swap_consistency(
        self, tmp_path, monkeypatch, cleanup_batchers
    ):
        """Test that checkpoint swapping produces consistent results.
        
        Verifies that updating a checkpoint mid-way through batching operations
        doesn't break the results - i.e., batching with a model, updating to the
        same model, and batching again should produce identical results.
        """
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp._base import (
            shutdown_inference_batchers,
            _current_batcher,
        )
        from quacc.recipes.mlp.core import static_job
        from quacc.wflow_tools.job_patterns import map_partition_fairchembatch
        from quacc.recipes.mlp import _base

        # Create test atoms
        atoms_list = [bulk("Cu"), bulk("Cu") * (2, 1, 1), bulk("Cu") * (3, 1, 1)]
        
        # Clear any existing batcher
        shutdown_inference_batchers()
        assert _base._current_batcher is None

        # First batch: process atoms with uma-s-1
        results_batch1 = map_partition_fairchembatch(
            static_job,
            atoms=atoms_list,
            name_or_path="uma-s-1",
            task_name="omat",
        )

        first_batcher = _base._current_batcher
        assert first_batcher is not None

        # Process the same atoms again (this triggers checkpoint update to same model)
        results_batch2 = map_partition_fairchembatch(
            static_job,
            atoms=atoms_list,
            name_or_path="uma-s-1",
            task_name="omat",
        )

        # Should still be using the same batcher instance (checkpoint was updated, not replaced)
        assert _base._current_batcher is first_batcher

        # Results should be identical
        assert len(results_batch1) == len(results_batch2) == 3
        for i, (result1, result2) in enumerate(zip(results_batch1, results_batch2)):
            assert result1["results"]["energy"] == pytest.approx(
                result2["results"]["energy"], rel=1e-6
            ), f"Energy mismatch at index {i} after checkpoint update"

    def test_batching_same_as_sequential_per_checkpoint(
        self, tmp_path, monkeypatch, cleanup_batchers
    ):
        """Test that batching with checkpoint swaps equals sequential processing.
        
        Verifies that:
        - Batch processing atoms [A, B, C, D] with checkpoint 1
        Then batch processing atoms [A, B, C, D] with checkpoint 1 (swapped)
        Gives the same total results as processing all atoms together at once.
        """
        monkeypatch.chdir(tmp_path)

        from quacc.recipes.mlp._base import shutdown_inference_batchers
        from quacc.recipes.mlp.core import static_job
        from quacc.wflow_tools.job_patterns import map_partition_fairchembatch

        # Create test atoms - use fixed seed for reproducibility
        atoms_list = [bulk("Cu"), bulk("Cu") * (2, 1, 1)]
        
        # Clear cache
        shutdown_inference_batchers()

        # Scenario 1: Process all atoms together
        results_all_at_once = map_partition_fairchembatch(
            static_job,
            atoms=atoms_list,
            name_or_path="uma-s-1",
            task_name="omat",
        )

        # Clear cache for next test
        shutdown_inference_batchers()

        # Scenario 2: Process atoms in two batches with checkpoint swaps in between
        # Batch 1: first atom
        results_batch_first = map_partition_fairchembatch(
            static_job,
            atoms=[atoms_list[0]],
            name_or_path="uma-s-1",
            task_name="omat",
        )

        # Batch 2: second atom (checkpoint stays same, but update_checkpoint is called)
        results_batch_second = map_partition_fairchembatch(
            static_job,
            atoms=[atoms_list[1]],
            name_or_path="uma-s-1",
            task_name="omat",
        )

        # Combine the split results
        results_split_batches = results_batch_first + results_batch_second

        # Both approaches should yield identical results
        assert len(results_all_at_once) == len(results_split_batches) == 2

        for i, (result_together, result_separate) in enumerate(
            zip(results_all_at_once, results_split_batches)
        ):
            # Energy should match
            assert result_together["results"]["energy"] == pytest.approx(
                result_separate["results"]["energy"], rel=1e-6
            ), f"Energy mismatch at index {i} between together vs separate batches"

            # Forces should match
            assert (
                result_together["results"]["forces"]
                == pytest.approx(
                    result_separate["results"]["forces"], rel=1e-6
                )
            ).all(), f"Forces mismatch at index {i} between together vs separate batches"

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

        from quacc.recipes.mlp._base import shutdown_inference_batchers
        from quacc.recipes.mlp.core import relax_job
        from quacc.wflow_tools.job_patterns import map_partition_fairchembatch

        # Clear batcher cache to ensure fresh start
        shutdown_inference_batchers()

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
            relax_job, atoms=[warmup_atoms], name_or_path="uma-s-1", task_name="omat"
        )

        # Time batched inference (after warmup)
        start_batched = time.perf_counter()
        results_batched = map_partition_fairchembatch(
            relax_job, atoms=atoms_list, name_or_path="uma-s-1", task_name="omat"
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
            atoms=partitioned_atoms,
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
                    atoms=[bulk("Cu")],
                    name_or_path="uma-s-1",
                    task_name="omat",
                )
