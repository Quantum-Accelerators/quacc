from __future__ import annotations

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


def test_merge_lists_preserves_case():
    defaults = ["%pal nprocs 4 end"]
    swaps = ['%moinp "/Data/Run1/Prev.gbw"', '%basis GTOName "MyBasis.BAS" end']
    assert merge_list_params(defaults, swaps) == [
        '%basis GTOName "MyBasis.BAS" end',
        '%moinp "/Data/Run1/Prev.gbw"',
        "%pal nprocs 4 end",
    ]

    # Duplicates and removals are still matched case-insensitively
    assert merge_list_params(
        ["B3LYP", "def2-SVP"], ["b3lyp", "#DEF2-svp", "TightSCF"]
    ) == ["b3lyp", "TightSCF"]
