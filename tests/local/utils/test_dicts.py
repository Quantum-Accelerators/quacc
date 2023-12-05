from quacc.utils.dicts import merge_dicts, remove_dict_nones,custom_deepcopy

def test_remove_dict_nones():
    d = {
        "output": {
            "output": {
                "test": [1, None],
                "test2": 1,
                "test3": None,
                "test4": {},
                "test5": [],
            }
        },
        "test": None,
    }
    d = remove_dict_nones(d)
    assert d == {
        "output": {"output": {"test": [1, None], "test2": 1, "test4": {}, "test5": []}}
    }

    # Additional test case 1: Nested dictionaries with None values
    d2 = {"a": None, "b": {"c": None, "d": 1}}
    d2 = remove_dict_nones(d2)
    assert d2 == {"b": {"d": 1}}

    # Additional test case 2: Empty dictionaries at different levels
    d3 = {"x": {}, "y": {"z": {}}}
    d3 = remove_dict_nones(d3)
    assert d3 == {}
    # Test custom_deepcopy with a list in remove_dict_nones context
    d4 = {"list": [1, [2, None], {"a": None}]}
    d4_copied = custom_deepcopy(d4)
    d4_copied = remove_dict_nones(d4_copied)
    assert d4_copied == {"list": [1, [2, None], {}]}

    # Test custom_deepcopy with a tuple in remove_dict_nones context
    d5 = {"tuple": (1, [2, None], {"a": None})}
    d5_copied = custom_deepcopy(d5)
    d5_copied = remove_dict_nones(d5_copied)
    assert d5_copied == {"tuple": (1, [2, None], {})}

def test_merge_dicts():
    defaults = {"a": 1, "b": {"a": 1, "b": 2}}
    calc_swaps = {"c": 3, "b": {"b": 3, "d": 1}}
    assert merge_dicts(defaults, calc_swaps) == {
        "a": 1,
        "b": {"a": 1, "b": 3, "d": 1},
        "c": 3,
    }

    # Additional test case 1: Merging dictionaries with missing keys
    merge_defaults = {"x": 10, "y": 20}
    merge_calc_swaps = {"y": 30, "z": 40}
    assert merge_dicts(merge_defaults, merge_calc_swaps) == {"x": 10, "y": 30, "z": 40}

def test_custom_deepcopy():
    # Test custom_deepcopy with a list containing None
    lst_with_none = [1, [2, None], {"a": None}]
    copied_list = custom_deepcopy(lst_with_none)
    assert copied_list == [1, [2, None], {}]

    # Test custom_deepcopy with a tuple containing None
    tpl_with_none = (1, [2, None], {"a": None})
    copied_tuple = custom_deepcopy(tpl_with_none)
    assert copied_tuple == (1, [2, None], {})

    # Test custom_deepcopy with a list containing nested tuples
    lst_with_tuples = [1, (2, [3, (4, None)]), {"a": ((5, [6, None]),)}]
    copied_lst_with_tuples = custom_deepcopy(lst_with_tuples)
    assert copied_lst_with_tuples == [1, (2, [3, (4, None)]), {"a": ((5, [6, None]),)}]

    # Test custom_deepcopy with a tuple containing nested lists
    tpl_with_lists = (1, [2, [3, [4, None]]], {"a": ([5, [6, None]],)})
    copied_tpl_with_lists = custom_deepcopy(tpl_with_lists)
    assert copied_tpl_with_lists == (1, [2, [3, [4, None]]], {"a": ([5, [6, None]],)})

    # Test custom_deepcopy with a list of tuples
    list_of_tuples = [(1, 2), (3, 4), (5, 6)]
    copied_list_of_tuples = custom_deepcopy(list_of_tuples)
    assert copied_list_of_tuples == [(1, 2), (3, 4), (5, 6)]

    # Test custom_deepcopy with a tuple of lists
    tuple_of_lists = ([1, 2], [3, 4], [5, 6])
    copied_tuple_of_lists = custom_deepcopy(tuple_of_lists)
    assert copied_tuple_of_lists == ([1, 2], [3, 4], [5, 6])