from quacc.util.dicts import merge_dicts, remove_dict_empties


def test_remove_dict_empties():
    d = {"output": {"output": {"test": [], "test2": 1}}}
    d = remove_dict_empties(d)
    assert d == {"output": {"output": {"test2": 1}}}

    d = {"output": {"output": {"test": {}, "test2": 1}}}
    d = remove_dict_empties(d)
    assert d == {"output": {"output": {"test2": 1}}}


def test_merge_dicts():
    defaults = {"a": 1, "b": {"a": 1, "b": 2}}
    calc_swaps = {"c": 3, "b": {"b": 3, "d": 1}}
    assert merge_dicts(defaults, calc_swaps) == {
        "a": 1,
        "b": {"a": 1, "b": 3, "d": 1},
        "c": 3,
    }
