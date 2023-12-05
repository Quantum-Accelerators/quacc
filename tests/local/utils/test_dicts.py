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
    # Additional test case 2: Merging dictionaries with nested dictionaries
    dict1 = {"a": {"b": 1}}
    dict2 = {"a": {"c": 2}}
    result = merge_dicts(dict1, dict2)
    assert result == {"a": {"b": 1, "c": 2}}

    # Additional test case 3: Merging dictionaries with nested lists
    dict3 = {"list1": [1, 2]}
    dict4 = {"list1": [3, 4]}
    result = merge_dicts(dict3, dict4)
    assert result == {"list1": [1, 2, 3, 4]}

    # Additional test case 4: Merging dictionaries with mixed types
    dict5 = {"mixed": {"a": [1, 2], "b": {"c": 3}}}
    dict6 = {"mixed": {"a": [4], "b": {"d": 5}}}
    result = merge_dicts(dict5, dict6)
    assert result == {"mixed": {"a": [1, 2, 4], "b": {"c": 3, "d": 5}}}

def test_custom_deepcopy_with_list_tuple():
    # Original list and tuple
    original_list = [1, [2, 3], {"a": 4}]
    original_tuple = (1, [2, 3], {"a": 4})

    # Copying the list and tuple
    copied_list = custom_deepcopy(original_list)
    copied_tuple = custom_deepcopy(original_tuple)

    # Modify the original list and tuple to check if the copies are indeed deep copies
    original_list[1].append(4)
    original_tuple[1].append(4)

    # Assertions to verify deep copying
    assert copied_list == [1, [2, 3], {"a": 4}]
    assert copied_tuple == (1, [2, 3], {"a": 4})
def test_merge_dicts_with_nested_dicts():
    # Dictionaries with a common key having dictionary values
    dict1 = {"common": {"key1": "value1"}, "unique1": "value2"}
    dict2 = {"common": {"key2": "value3"}, "unique2": "value4"}

    # Merging the dictionaries
    merged_dict = merge_dicts(dict1, dict2)

    # Expected result: the nested dictionaries under "common" should be merged
    assert merged_dict == {
        "common": {"key1": "value1", "key2": "value3"},
        "unique1": "value2",
        "unique2": "value4"
    }
