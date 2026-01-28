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
        "orb": -4.093477725982666,
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
        "orb": -32.7361946105957,
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
        "orb": -32.73428726196289,
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
def test_formation_energy_error_without_omat(tmp_path, monkeypatch):
    """Test that formation energy raises error when not using omat task."""
    monkeypatch.chdir(tmp_path)
    _set_dtype(32)

    from huggingface_hub.utils._auth import get_token

    if not get_token():
        pytest.skip("HuggingFace token not available for FAIRChem")

    # Try to use formation energy without omat task - should raise
    atoms = bulk("Cu")
    with pytest.raises(ValueError, match="task_name='omat'"):
        static_job(
            atoms,
            method="fairchem",
            name_or_path="uma-s-1",
            task_name="omol",  # Wrong task
            use_formation_energy=True,
        )


@pytest.mark.skipif(find_spec("fairchem") is None, reason="fairchem not installed")
def test_formation_energy_error_without_fairchem(tmp_path, monkeypatch):
    """Test that formation energy raises error when not using fairchem."""
    monkeypatch.chdir(tmp_path)

    if "mace-mp" not in methods:
        pytest.skip("mace-mp not available")

    _set_dtype(64)

    # Try to use formation energy with non-FAIRChem method - should raise
    atoms = bulk("Cu")
    with pytest.raises(ValueError, match="FAIRChem UMA"):
        static_job(atoms, method="mace-mp", use_formation_energy=True)


@pytest.mark.skipif(find_spec("fairchem") is None, reason="fairchem not installed")
def test_static_job_formation_energy_with_kwargs(tmp_path, monkeypatch):
    """Test that formation_energy_kwargs are properly passed through."""
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
        formation_energy_kwargs={},  # Pass empty dict as formation_energy_kwargs
    )
    # Should succeed and compute formation energy
    assert abs(output["results"]["energy"]) < 0.1
    assert np.shape(output["results"]["forces"]) == (1, 3)
