from quacc.util.json import jsanitize, jdesanitize
from ase.build import bulk
from ase.io.jsonio import encode
from pymatgen.io.ase import AseAtomsAdaptor
import numpy as np


def test_json():
    atoms = bulk("Cu")
    assert jsanitize(atoms) == encode(atoms)
    assert jdesanitize(jsanitize(atoms)) == atoms
    struct = AseAtomsAdaptor.get_structure(bulk("Mg"))
    assert jdesanitize(jsanitize(struct)) == struct

    atoms2 = bulk("Al")
    atoms2.info = {"foo": "bar"}
    assert jdesanitize(jsanitize(atoms2)).info == {"foo": "bar"}

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
    assert jdesanitize(jsanitize(test_set["test"])) == [1, 2, 3]
    assert jdesanitize(jsanitize(test_set["test2"])) == atoms2
    assert jdesanitize(jsanitize(test_set["test3"])) == struct
    assert jdesanitize(jsanitize(test_set["test4"])) == "hi"
    assert jdesanitize(jsanitize(test_set["test5"])) == [1, 2, 3]
    assert jdesanitize(jsanitize(test_set["test6"])) == 1
    assert jdesanitize(jsanitize(test_set["test7"])) == [[1, 2, 3], [4, 5, 6]]
    assert jdesanitize(jsanitize(test_set["test8"])) == True
    assert jdesanitize(jsanitize(test_set["test9"])) == None
