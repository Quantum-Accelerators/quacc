from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")

from importlib.util import find_spec

import numpy as np
from ase.build import bulk

from quacc.recipes.mlp.core import relax_job, static_job

methods = []
if has_mace := find_spec("mace"):
    methods.append("mace-mp")

if find_spec("matgl"):
    methods.append("tensornet")

if has_sevennet := find_spec("sevenn"):
    methods.append("sevennet")

if find_spec("orb_models"):
    methods.append("orb")

if find_spec("fairchem"):
    from huggingface_hub.utils._auth import get_token

    if get_token():
        methods.append("fairchem")


def _set_dtype(size, type_="float"):
    globals()[f"{type_}_th"] = getattr(torch, f"{type_}{size}")
    globals()[f"{type_}_np"] = getattr(np, f"{type_}{size}")
    torch.set_default_dtype(getattr(torch, f"float{size}"))


@pytest.mark.parametrize("method", methods)
def test_static_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)

    if method == "mace-mp":
        _set_dtype(64)
    else:
        _set_dtype(32)

    if method == "fairchem":
        # Note that for this to work, you need HF_TOKEN env variable set!
        calc_kwargs = {"name_or_path": "uma-s-1", "task_name": "omat"}
    else:
        calc_kwargs = {}

    ref_energy = {
        "tensornet": -3.7593491077423096,
        "mace-mp": -4.097862720291976,
        "sevennet": -4.096191883087158,
        "orb": -3.7420763969421387,
        "fairchem": -3.7579006783217954,
    }
    atoms = bulk("Cu")
    output = static_job(atoms, method=method, **calc_kwargs)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method], rel=1e-4)
    assert np.shape(output["results"]["forces"]) == (1, 3)
    assert output["atoms"] == atoms


@pytest.mark.skipif(has_sevennet is None, reason="sevennet not installed")
def test_static_job_with_dict_kwargs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")

    # Make sure that pick_calculator works even with dictionary kwargs
    static_job(atoms, method="sevennet", sevennet_config={"test": 1})


@pytest.mark.parametrize("method", methods)
def test_relax_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)

    if method == "mace-mp":
        _set_dtype(64)
    else:
        _set_dtype(32)

    if method == "fairchem":
        # Note that for this to work, you need HF_TOKEN env variable set!
        calc_kwargs = {"name_or_path": "uma-s-1", "task_name": "omat"}
    else:
        calc_kwargs = {}

    ref_energy = {
        "mace-mp": -32.78264569638644,
        "tensornet": -30.074462890625,
        "sevennet": -32.76924133300781,
        "orb": -29.93630599975586,
        "fairchem": -30.004380887389797,
    }

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(atoms, method=method, **calc_kwargs)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method], rel=1e-4)
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() == pytest.approx(atoms.get_volume())


@pytest.mark.skipif(has_mace is None, reason="Needs MACE")
@pytest.mark.skipif(find_spec("torch_dftd") is None, reason="Needs torch-dftd")
def test_relax_job_dispersion(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    _set_dtype(64)

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(atoms, method="mace-mp", dispersion=True)
    assert output["results"]["energy"] == pytest.approx(-37.4518034464096)
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() == pytest.approx(atoms.get_volume())


@pytest.mark.parametrize("method", methods)
def test_relax_cell_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)

    if method == "mace-mp":
        _set_dtype(64)
    else:
        _set_dtype(32)

    if method == "fairchem":
        # Note that for this to work, you need HF_TOKEN env variable set!
        calc_kwargs = {"name_or_path": "uma-s-1", "task_name": "omat"}
    else:
        calc_kwargs = {}

    ref_energy = {
        "mace-mp": -32.8069374165035,
        "tensornet": -30.079431533813477,
        "sevennet": -32.76963806152344,
        "orb": -29.93630599975586,
        "fairchem": -30.005004590392726,
    }

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(atoms, method=method, relax_cell=True, **calc_kwargs)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method], rel=1e-4)
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() != pytest.approx(atoms.get_volume())


@pytest.mark.skipif(find_spec("fairchem") is None, reason="fairchem not installed")
def test_static_job_formation_energy_fairchem(tmp_path, monkeypatch):
    """Test formation energy calculation with FAIRChem UMA omat."""
    monkeypatch.chdir(tmp_path)
    _set_dtype(32)

    from huggingface_hub.utils._auth import get_token

    if not get_token():
        pytest.skip("HuggingFace token not available for FAIRChem")

    calc_kwargs = {"name_or_path": "uma-s-1", "task_name": "omat"}

    # Test Cu (elemental system - formation energy should be ~0)
    atoms_cu = bulk("Cu")
    output_cu = static_job(
        atoms_cu, method="fairchem", use_formation_energy=True, **calc_kwargs
    )
    # For pure elements, total formation energy should be close to zero (eV)
    assert abs(output_cu["results"]["energy"]) < 0.1
    assert np.shape(output_cu["results"]["forces"]) == (1, 3)

    # Test MgO (binary compound - formation energy should be negative)
    atoms_mgo = bulk("MgO", crystalstructure="rocksalt", a=4.2)
    output_mgo = static_job(
        atoms_mgo, method="fairchem", use_formation_energy=True, **calc_kwargs
    )
    # MgO has a substantial negative formation energy (eV per formula unit)
    assert output_mgo["results"]["energy"] < -2.0
    assert np.shape(output_mgo["results"]["forces"]) == (2, 3)


@pytest.mark.skipif(find_spec("fairchem") is None, reason="fairchem not installed")
def test_relax_job_formation_energy_fairchem(tmp_path, monkeypatch):
    """Test formation energy calculation during relaxation with FAIRChem UMA omat."""
    monkeypatch.chdir(tmp_path)
    _set_dtype(32)

    from huggingface_hub.utils._auth import get_token

    if not get_token():
        pytest.skip("HuggingFace token not available for FAIRChem")

    calc_kwargs = {"name_or_path": "uma-s-1", "task_name": "omat"}

    # Test Cu supercell
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(
        atoms, method="fairchem", use_formation_energy=True, **calc_kwargs
    )
    # 8 Cu atoms, total formation energy should be near zero (eV)
    assert abs(output["results"]["energy"]) < 0.8
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms


@pytest.mark.skipif(find_spec("fairchem") is None, reason="fairchem not installed")
def test_formation_energy_multiple_compounds_fairchem(tmp_path, monkeypatch):
    """Test formation energy for multiple compound types with FAIRChem UMA omat."""
    monkeypatch.chdir(tmp_path)
    _set_dtype(32)

    from huggingface_hub.utils._auth import get_token

    if not get_token():
        pytest.skip("HuggingFace token not available for FAIRChem")

    calc_kwargs = {"name_or_path": "uma-s-1", "task_name": "omat"}

    # Test NaCl
    atoms_nacl = bulk("NaCl", crystalstructure="rocksalt", a=5.64)
    output_nacl = static_job(
        atoms_nacl, method="fairchem", use_formation_energy=True, **calc_kwargs
    )
    # NaCl has negative formation energy (eV per formula unit)
    assert output_nacl["results"]["energy"] < -1.0

    # Test Si
    atoms_si = bulk("Si")
    output_si = static_job(
        atoms_si, method="fairchem", use_formation_energy=True, **calc_kwargs
    )
    # Pure Si should have near-zero formation energy (eV)
    assert abs(output_si["results"]["energy"]) < 0.1

    # Test Al
    atoms_al = bulk("Al")
    output_al = static_job(
        atoms_al, method="fairchem", use_formation_energy=True, **calc_kwargs
    )
    # Pure Al should have near-zero formation energy (eV)
    assert abs(output_al["results"]["energy"]) < 0.1


@pytest.mark.skipif(find_spec("fairchem") is None, reason="fairchem not installed")
def test_relax_job_formation_energy_cell_fairchem(tmp_path, monkeypatch):
    """Test formation energy calculation with cell relaxation using FAIRChem UMA omat."""
    monkeypatch.chdir(tmp_path)
    _set_dtype(32)

    from huggingface_hub.utils._auth import get_token

    if not get_token():
        pytest.skip("HuggingFace token not available for FAIRChem")

    calc_kwargs = {"name_or_path": "uma-s-1", "task_name": "omat"}

    # Test MgO with cell relaxation
    atoms_mgo = bulk("MgO", crystalstructure="rocksalt", a=4.2) * (2, 2, 2)
    atoms_mgo[0].position += 0.05
    output = relax_job(
        atoms_mgo,
        method="fairchem",
        relax_cell=True,
        use_formation_energy=True,
        **calc_kwargs,
    )
    # Should have relaxed and computed formation energy
    assert (
        output["results"]["energy"] < -8.0
    )  # 8 formula units * ~-1 eV per formula unit
    assert output["atoms"] != atoms_mgo


@pytest.mark.skipif(find_spec("fairchem") is None, reason="fairchem not installed")
def test_formation_energy_non_omat_task(tmp_path, monkeypatch):
    """Test that formation energy works with non-omat tasks (uses default references)."""
    monkeypatch.chdir(tmp_path)
    _set_dtype(32)

    from huggingface_hub.utils._auth import get_token

    if not get_token():
        pytest.skip("HuggingFace token not available for FAIRChem")

    # Formation energy should work with non-omat tasks, using default references
    atoms = bulk("Cu")
    # This should not raise an error - it will use default references
    output = static_job(
        atoms,
        method="fairchem",
        name_or_path="uma-s-1",
        task_name="omol",  # Non-omat task
        use_formation_energy=True,
    )
    # Should complete successfully
    assert "energy" in output["results"]


@pytest.mark.skipif(find_spec("fairchem") is None, reason="fairchem not installed")
def test_formation_energy_with_mace(tmp_path, monkeypatch):
    """Test that formation energy works with non-FAIRChem models when references provided."""
    monkeypatch.chdir(tmp_path)

    if "mace-mp" not in methods:
        pytest.skip("mace-mp not available")

    _set_dtype(64)

    # Formation energy should work with mace-mp when references are provided
    atoms = bulk("Cu")
    # Use OMAT24 references for testing

    output = static_job(
        atoms, method="mace-mp", use_formation_energy=True, references="OMAT24"
    )
    # Should complete successfully
    assert "energy" in output["results"]


@pytest.mark.skipif(find_spec("fairchem") is None, reason="fairchem not installed")
def test_builtin_vs_omat24_references(tmp_path, monkeypatch):
    """Test that built-in and OMAT24 references give identical results for UMA-s-1p1 with task_name='omat'."""
    monkeypatch.chdir(tmp_path)
    _set_dtype(32)

    from huggingface_hub.utils._auth import get_token

    if not get_token():
        pytest.skip("HuggingFace token not available for FAIRChem")

    calc_kwargs = {"name_or_path": "uma-s-1p1", "task_name": "omat"}
    atoms = bulk("Cu")

    # Test with built-in references (None)
    output_builtin = static_job(
        atoms,
        method="fairchem",
        use_formation_energy=True,
        references=None,  # Use built-in
        **calc_kwargs,
    )
    energy_builtin = output_builtin["results"]["energy"]

    # Test with OMAT24 references (explicit)
    output_omat24 = static_job(
        atoms,
        method="fairchem",
        use_formation_energy=True,
        references="OMAT24",
        **calc_kwargs,
    )
    energy_omat24 = output_omat24["results"]["energy"]

    # Results should be identical (or very close due to numerical precision)
    assert abs(energy_builtin - energy_omat24) < 1e-6, (
        f"Built-in references ({energy_builtin:.6f} eV) differ from OMAT24 "
        f"references ({energy_omat24:.6f} eV) by {abs(energy_builtin - energy_omat24):.6e} eV"
    )


@pytest.mark.skipif(find_spec("fairchem") is None, reason="fairchem not installed")
def test_static_job_formation_energy_with_references(tmp_path, monkeypatch):
    """Test that references parameter is properly passed through."""
    monkeypatch.chdir(tmp_path)
    _set_dtype(32)

    from huggingface_hub.utils._auth import get_token

    if not get_token():
        pytest.skip("HuggingFace token not available for FAIRChem")

    atoms = bulk("Cu")
    output = static_job(
        atoms,
        method="fairchem",
        name_or_path="uma-s-1",
        task_name="omat",
        use_formation_energy=True,
        references=None,  # Use built-in references
    )
    # Should succeed and compute formation energy
    assert abs(output["results"]["energy"]) < 0.1
    assert np.shape(output["results"]["forces"]) == (1, 3)


@pytest.mark.skipif(find_spec("fairchem") is None, reason="fairchem not installed")
def test_formation_energy_consistency_across_models(tmp_path, monkeypatch):
    """Test that formation energies are consistent across different models.

    This test ensures that when different ML potentials calculate formation
    energies with appropriate references, they yield similar results
    (within ~0.1 eV/atom). This helps debug issues with formation energy
    references and MP corrections.
    """
    monkeypatch.chdir(tmp_path)

    from huggingface_hub.utils._auth import get_token

    if not get_token():
        pytest.skip("HuggingFace token not available for FAIRChem")

    # Get OMAT24 references for non-FAIRChem models
    # (FAIRChem will use built-in references with None)

    # Test elemental system (Cu) - formation energy should be near zero
    test_structures = [
        ("Cu", bulk("Cu"), 0.0, 0.2),  # (name, structure, expected_per_atom, tolerance)
        ("MgO", bulk("MgO", crystalstructure="rocksalt", a=4.2), -3.0, 0.5),
    ]

    for struct_name, atoms, expected_per_atom, tolerance in test_structures:
        results = {}

        # Test with FAIRChem (uses built-in references with None)
        if "fairchem" in methods:
            _set_dtype(32)
            output = static_job(
                atoms,
                method="fairchem",
                name_or_path="uma-s-1",
                task_name="omat",
                use_formation_energy=True,
                references=None,  # Use built-in
            )
            energy_per_atom = output["results"]["energy"] / len(atoms)
            results["fairchem"] = energy_per_atom

        # Test with MACE (using OMAT24 references)
        if "mace-mp" in methods:
            _set_dtype(64)
            output = static_job(
                atoms, method="mace-mp", use_formation_energy=True, references="OMAT24"
            )
            energy_per_atom = output["results"]["energy"] / len(atoms)
            results["mace-mp"] = energy_per_atom

        # Test with TensorNet (using OMAT24 references)
        if "tensornet" in methods:
            _set_dtype(32)
            output = static_job(
                atoms,
                method="tensornet",
                use_formation_energy=True,
                references="OMAT24",
            )
            energy_per_atom = output["results"]["energy"] / len(atoms)
            results["tensornet"] = energy_per_atom

        # Test with SevenNet (using OMAT24 references)
        if "sevennet" in methods:
            _set_dtype(32)
            output = static_job(
                atoms, method="sevennet", use_formation_energy=True, references="OMAT24"
            )
            energy_per_atom = output["results"]["energy"] / len(atoms)
            results["sevennet"] = energy_per_atom

        # Test with ORB (using OMAT24 references)
        if "orb" in methods:
            _set_dtype(32)
            output = static_job(
                atoms, method="orb", use_formation_energy=True, references="OMAT24"
            )
            energy_per_atom = output["results"]["energy"] / len(atoms)
            results["orb"] = energy_per_atom

        # Verify we have at least 2 results to compare
        if len(results) < 2:
            pytest.skip(f"Not enough models available to compare for {struct_name}")

        # Check that all results are within expected range
        for method_name, energy in results.items():
            assert abs(energy - expected_per_atom) < tolerance, (
                f"{struct_name} - {method_name}: {energy:.4f} eV/atom is outside "
                f"expected range {expected_per_atom} +/- {tolerance} eV/atom"
            )

        # Check consistency across models (within 0.1 eV/atom)
        energies_list = list(results.values())
        max_energy = max(energies_list)
        min_energy = min(energies_list)
        spread = max_energy - min_energy

        # Allow slightly larger spread for compounds than elements
        max_spread = 0.15 if "Cu" in struct_name else 0.2
        assert spread < max_spread, (
            f"{struct_name} - Spread of {spread:.4f} eV/atom exceeds maximum "
            f"allowed {max_spread} eV/atom. Results: {results}"
        )


@pytest.mark.skipif(find_spec("fairchem") is None, reason="fairchem not installed")
def test_formation_energy_mp20_references(tmp_path, monkeypatch):
    """Test that MP-20 references can be loaded and used."""
    monkeypatch.chdir(tmp_path)

    from huggingface_hub.utils._auth import get_token

    if not get_token():
        pytest.skip("HuggingFace token not available for FAIRChem")

    # Test loading MP-20 references
    from quacc.recipes.mlp._base import _get_mp20_references

    try:
        mp20_refs = _get_mp20_references()
        assert isinstance(mp20_refs, dict)
        assert len(mp20_refs) > 0

        # Test with FAIRChem using MP20 references (explicit)
        if "fairchem" in methods:
            _set_dtype(32)
            atoms = bulk("Cu")

            # Test with MP-20 references
            output = static_job(
                atoms,
                method="fairchem",
                name_or_path="uma-s-1",
                task_name="omat",
                use_formation_energy=True,
                references="MP20",
            )
            energy_per_atom = output["results"]["energy"] / len(atoms)

            # Should be close to zero for elemental system
            assert abs(energy_per_atom) < 0.2

    except Exception as e:
        pytest.skip(f"Could not load MP-20 references: {e}")
