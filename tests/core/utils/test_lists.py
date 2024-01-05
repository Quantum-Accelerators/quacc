from quacc.utils.lists import merge_list_params


def test_merge_lists():
    a = ["a", "b", "c"]
    b = ["b", "c", "d"]
    c = ["c", "d", "e"]
    assert merge_list_params(a, b, c) == ["a", "b", "c", "d", "e"]

    a = ["a", "b", "c"]
    b = ["b", "c", "d"]
    c = ["C", "d", "e"]
    assert merge_list_params(a, b, c, case_insensitive=False) == [
        "C",
        "a",
        "b",
        "c",
        "d",
        "e",
    ]

    a = ["a", "b", "c"]
    b = ["b", "c", "d"]
    c = ["#c", "d", "e"]
    assert merge_list_params(a, b, c) == ["a", "b", "d", "e"]

    a = ["a", "b", "c"]
    b = ["b", "c", "d"]
    c = ["#C", "d", "e"]
    assert merge_list_params(a, b, c, case_insensitive=False) == [
        "a",
        "b",
        "c",
        "d",
        "e",
    ]

    a = ["wb97x-d3bj", "def2-tzvp", "opt", "slowconv", "normalprint", "xyzfile"]
    b = ["def2-svp", "#def2-tzvp", "#Wb97x-d3bj", "hf"]
    assert merge_list_params(a, b) == [
        "def2-svp",
        "hf",
        "normalprint",
        "opt",
        "slowconv",
        "xyzfile",
    ]
