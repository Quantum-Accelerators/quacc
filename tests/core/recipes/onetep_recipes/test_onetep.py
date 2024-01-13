from pathlib import Path

from ase.build import molecule

from quacc import SETTINGS
from quacc.recipes.onetep.core import static_job
from quacc.utils.files import copy_decompress_files

DEFAULT_SETTINGS = SETTINGS.model_copy()
DATA_DIR = Path(__file__).parent / "data"


def test_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "H.usp.gz"], tmp_path)

    keywords = {"pseudo_path": str(tmp_path)}
    pseudopotentials = {"H": "H.usp"}

    atoms = molecule("H2")
    atoms.set_cell([10, 10, 10])
    atoms.center()

    output = static_job(
        atoms,
        keywords=keywords,
        pseudopotentials=pseudopotentials,
        ngwf_radius=6.0,
        ngwf_count=1,
    )

    assert output["natoms"] == len(atoms)
    assert output["atoms"] == atoms
    assert output["parameters"]["keywords"] == {
        "output_detail": "verbose",
        "do_properties": True,
        "cutoff_energy": "600 eV",
        "task": "SinglePoint",
        "pseudo_path": str(tmp_path),
    }
