from ase.build import molecule
from quacc import SETTINGS
from pathlib import Path
from quacc.recipes.gaussian.core import relax_job, static_job
import pytest
from quacc.utils.files import copy_decompress_files

DEFAULT_SETTINGS = SETTINGS.model_copy()
DATA_DIR = Path(__file__).parent / "data"


def test_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    pseudopotentials={"H":"H_NCP19_PBE_OTF.usp"}
    pp_dir = Path(__file__).parent
    keywords = {
            'output_detail': 'verbose',
            'do_properties': True,
            'cutoff_energy' : '300 eV',
            'task' : 'SinglePoint'

}
    copy_decompress_files([pp_dir / "H_NCP19_PBE_OTF.usp"], tmp_path)
    atoms = molecule("H2")

    output = static_job(atoms, pseudo_path=pseudo_path, pseudopotentials = pseudopotentials)
    assert output["natoms"] == len(atoms)
    assert output["elements"] == "H"
    assert output["parameters"]["keywords"] == keywords