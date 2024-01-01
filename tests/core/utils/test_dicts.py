from quacc import Remove
from quacc.utils.dicts import recursive_dict_merge, remove_dict_entries


def test_remove_dict_entries():
    d = {
        "output": {
            "output": {
                "test": [1, Remove],
                "test2": 1,
                "test3": Remove,
                "test4": {},
                "test5": [],
            }
        },
        "test": Remove,
    }
    d = remove_dict_entries(d)
    assert d == {
        "output": {
            "output": {"test": [1, Remove], "test2": 1, "test4": {}, "test5": []}
        }
    }


def test_recursive_dict_merge():
    defaults = {"a": 1, "b": {"a": 1, "b": 2}}
    calc_swaps = {"c": 3, "b": {"b": 3, "d": 1}}
    assert recursive_dict_merge(defaults, calc_swaps) == {
        "a": 1,
        "b": {"a": 1, "b": 3, "d": 1},
        "c": 3,
    }
