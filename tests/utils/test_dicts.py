from quacc.utils.dicts import remove_dict_empties


def test_remove_dict_empties():
    d = {"output": {"output": {"test": [], "test2": 1}}}
    d = remove_dict_empties(d)
    assert d == {"output": {"output": {"test2": 1}}}

    d = {"output": {"output": {"test": {}, "test2": 1}}}
    d = remove_dict_empties(d)
    assert d == {"output": {"output": {"test2": 1}}}
