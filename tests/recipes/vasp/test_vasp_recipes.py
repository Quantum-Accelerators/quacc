import os
from pathlib import Path
from shutil import copy

from ase.build import bulk
from jobflow.managers.local import run_locally

from quacc.recipes.vasp.core import RelaxMaker, StaticMaker
from quacc.recipes.vasp.slabs import (
    BulkToSlabMaker,
    SlabRelaxMaker,
    SlabStaticMaker,
    SlabToAdsSlabMaker,
)
from quacc.util.json import jsonify

FILE_DIR = Path(__file__).resolve().parent
run1 = os.path.join(FILE_DIR, "..", "..", "schemas", "vasp_run1")

# Copy VASP files to current directory so the "output" can be parsed
def setup_module():
    for f in os.listdir(run1):
        copy(os.path.join(run1, f), os.path.join(FILE_DIR, f))


# Remove copied VASP files after tests are done
def teardown_module():
    for f in os.listdir(run1):
        os.remove(os.path.join(FILE_DIR, f))


def test_core():
    atoms_json = jsonify(bulk("Cu"))

    job = StaticMaker().make(atoms_json)
    output = run_locally(job)
    assert output is not None

    job = RelaxMaker().make(atoms_json)
    output = run_locally(job)
    assert output is not None


def test_slabs():
    atoms_json = jsonify(bulk("Cu"))

    job = SlabStaticMaker().make(atoms_json)
    output = run_locally(job)
    assert output is not None

    job = SlabRelaxMaker().make(atoms_json)
    output = run_locally(job)
    assert output is not None

    job = BulkToSlabMaker().make(atoms_json)
    output = run_locally(job)
    assert output is not None

    atoms = bulk("Cu") * (2, 2, 2)
    atoms.center(vacuum=10, axis=2)
    slab_json = jsonify(atoms)
    job = SlabToAdsSlabMaker().make(slab_json)
    output = run_locally(job)
    assert output is not None
