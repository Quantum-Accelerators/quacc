from ase.build import bulk
from quacc.recipes.vasp.core import RelaxMaker, StaticMaker
from quacc.util.json import jsonify
from jobflow.managers.local import run_locally
from pathlib import Path
from shutil import copy
import os

FILE_DIR = Path(__file__).resolve().parent
run1 = os.path.join(FILE_DIR, "..", "..", "schemas", "vasp_run1")


def setup_module():
    for f in os.listdir(run1):
        copy(os.path.join(run1, f), os.path.join(FILE_DIR, f))


def teardown_module():
    for f in os.listdir(run1):
        os.remove(os.path.join(FILE_DIR, f))


def test_static_maker():

    atoms_json = jsonify(bulk("Cu"))
    job = StaticMaker().make(atoms_json)
    output = run_locally(job)
    assert output is not None