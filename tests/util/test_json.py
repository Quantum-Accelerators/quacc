from quacc.util.json import jsonify, unjsonify
from ase.build import bulk
from ase.io.jsonio import encode
from pymatgen.io.ase import AseAtomsAdaptor
import numpy as np


def test_json():
    atoms = bulk("Cu")
    assert jsonify(atoms) == encode(atoms)
    assert unjsonify(jsonify(atoms)) == atoms
    struct = AseAtomsAdaptor.get_structure(bulk("Mg"))
    assert unjsonify(jsonify(struct)) == struct

    atoms2 = bulk("Al")
    atoms2.info = {"foo": "bar"}
    assert unjsonify(jsonify(atoms2)).info == {"foo": "bar"}

    test_set = {
        "test": np.array([1, 2, 3]),
        "test2": atoms2,
        "test3": struct,
        "test4": "hi",
        "test5": [1, 2, 3],
        "test6": 1,
        "test7": [[1, 2, 3], np.array([4, 5, 6])],
        "test8": True,
        "test9": None,
    }
    assert unjsonify(jsonify(test_set["test"])) == [1, 2, 3]
    assert unjsonify(jsonify(test_set["test2"])) == atoms2
    assert unjsonify(jsonify(test_set["test3"])) == struct
    assert unjsonify(jsonify(test_set["test4"])) == "hi"
    assert unjsonify(jsonify(test_set["test5"])) == [1, 2, 3]
    assert unjsonify(jsonify(test_set["test6"])) == 1
    assert unjsonify(jsonify(test_set["test7"])) == [[1, 2, 3], [4, 5, 6]]
    assert unjsonify(jsonify(test_set["test8"])) == True
    assert unjsonify(jsonify(test_set["test9"])) == None
