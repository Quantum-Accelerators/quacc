from __future__ import annotations

import numpy as np
import pytest
from ase.build import bulk

from quacc.recipes.mlp.adsorbml import bulk_to_surfaces_to_adsorbml
from quacc.recipes.mlp.core import relax_job

torch = pytest.importorskip("torch")

from importlib.util import find_spec

methods = []
if has_mace := find_spec("mace"):
    methods.append("mace-mp-0")

if has_matgl := find_spec("matgl"):
    methods.append("m3gnet")

if has_chgnet := find_spec("chgnet"):
    methods.append("chgnet")

if has_sevennet := find_spec("sevenn"):
    methods.append("sevennet")

if has_orb := find_spec("orb_models"):
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
def test_total_energy_adsorbml(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)

    if method == "mace-mp-0":
        _set_dtype(64)
    else:
        _set_dtype(32)

    if method == "fairchem":
        calc_kwargs = {
            "model_name": "EquiformerV2-31M-OMAT24-MP-sAlex",
            "local_cache": "./fairchem_checkpoint_cache/",
            "seed": 42,
        }
    else:
        calc_kwargs = {}

    ref_CO_Cu111_adsorption_energy = {
        "chgnet": -1.5,
        "m3gnet": -0.5,
        "mace-mp-0": -0.5,
        "sevennet": -0.98,
        "orb": -0.5,
        "fairchem": -0.5,
    }

    bulk_atoms = bulk("Cu")

    outputs = bulk_to_surfaces_to_adsorbml(
        bulk_atoms,
        adsorbates_kwargs=[{"adsorbate_smiles_from_db": "*CO"}],
        multiple_adsorbate_slab_config_kwargs={"num_configurations": 1},
        ml_relax_job=relax_job,
        slab_validate_job=relax_job,
        adslab_validate_job=relax_job,
        gas_validate_job=relax_job,
        bulk_relax_job=relax_job,
        job_params={
            "all": dict(method=method, **calc_kwargs, opt_params={"max_steps": 500})
        },
        max_miller=1,
        num_to_validate_with_DFT=1,
        reference_ml_energies_to_gas_phase=True,
        relax_bulk=True,
    )
    index_111_facet = np.where(
        [result["slab"]["slab_metadata"]["millers"] == (1, 1, 1) for result in outputs]
    )[0][0]

    assert len(outputs) > 0
    assert outputs[index_111_facet]["adslabs"][0]["results"][
        "referenced_adsorption_energy"
    ]["adsorption_energy"] == pytest.approx(
        ref_CO_Cu111_adsorption_energy[method], abs=0.1
    )


@pytest.mark.parametrize("fairchem_checkpoint", ["EquiformerV2-31M-S2EF-OC20-All+MD"])
def test_referenced_energy_mlp(tmp_path, monkeypatch, fairchem_checkpoint):
    monkeypatch.chdir(tmp_path)

    calc_kwargs = {
        "model_name": fairchem_checkpoint,
        "local_cache": "./fairchem_checkpoint_cache/",
        "seed": 42,
    }

    ref_CO_Cu111_adsorption_energy = {"EquiformerV2-31M-S2EF-OC20-All+MD": -0.5}

    bulk_atoms = bulk("Cu")

    outputs = bulk_to_surfaces_to_adsorbml(
        bulk_atoms,
        adsorbates_kwargs=[{"adsorbate_smiles_from_db": "*CO"}],
        multiple_adsorbate_slab_config_kwargs={"num_configurations": 1},
        ml_relax_job=relax_job,
        slab_validate_job=relax_job,
        adslab_validate_job=relax_job,
        gas_validate_job=relax_job,
        bulk_relax_job=relax_job,
        job_params={"all": dict(method="fairchem", **calc_kwargs)},
        max_miller=1,
        num_to_validate_with_DFT=1,
        reference_ml_energies_to_gas_phase=False,
        relax_bulk=False,
    )
    index_111_facet = np.where(
        [result["slab"]["slab_metadata"]["millers"] == (1, 1, 1) for result in outputs]
    )[0][0]

    assert len(outputs) > 0
    assert outputs[index_111_facet]["adslabs"][0]["results"]["energy"] == pytest.approx(
        ref_CO_Cu111_adsorption_energy[fairchem_checkpoint], abs=0.1
    )
