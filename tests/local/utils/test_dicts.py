from quacc.utils.dicts import merge_dicts, remove_dict_nones,custom_deepcopy

def test_remove_dict_nones():
    d = {
        "output": {
            "output": {
                "test": [1, None],
                "test2": 1,
                "test3": None,
            }
        },
        "test": None,
    }
    d = remove_dict_nones(d)
    assert d == {
        "output": {"output": {"test": [1, None], "test2": 1}}
    }

    # Additional test case 1: Nested dictionaries with None values
    d2 = {"a": None, "b": {"c": None, "d": 1}}
    d2 = remove_dict_nones(d2)
    assert d2 == {"b": {"d": 1}}

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

    # Additional test case 2: Merging dictionaries with nested dictionaries
    dict1 = {"a": {"b": 1}}
    dict2 = {"a": {"c": 2}}
    result = merge_dicts(dict1, dict2)
    assert result == {"a": {"b": 1, "c": 2}}

def test_custom_deepcopy_with_simple_list_tuple():
    # Simple list and tuple
    original_list = [1, 2, 3]
    original_tuple = (4, 5, 6)

    # Copy them using custom_deepcopy
    copied_list = custom_deepcopy(original_list)
    copied_tuple = custom_deepcopy(original_tuple)

    # Assertions
    assert copied_list == original_list
    assert copied_tuple == original_tuple
    # Ensure they are not the same objects (deep copy check) for the list only
    assert copied_list is not original_list

def test_custom_deepcopy_with_simple_list_tuple():
    # Simple list and tuple
    original_list = [1, 2, 3]
    original_tuple = (4, 5, 6)

    # Copy them using custom_deepcopy
    copied_list = custom_deepcopy(original_list)
    copied_tuple = custom_deepcopy(original_tuple)

    # Assertions
    assert copied_list == original_list
    assert copied_tuple == original_tuple
    # Ensure they are not the same objects (deep copy check)
    assert copied_list is not original_list
    assert copied_tuple is not original_tuple
