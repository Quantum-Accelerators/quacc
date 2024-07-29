from __future__ import annotations

import logging

import pytest

from quacc import Remove
from quacc.utils.dicts import finalize_dict, recursive_dict_merge, remove_dict_entries

LOGGER = logging.getLogger(__name__)
LOGGER.propagate = True


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
    d = remove_dict_entries(d, Remove)
    assert d == {
        "output": {
            "output": {"test": [1, Remove], "test2": 1, "test4": {}, "test5": []}
        }
    }


def test_remove_dict_entries2():
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
    d = remove_dict_entries(d, None)
    assert d == {
        "output": {"output": {"test": [1, None], "test2": 1, "test4": {}, "test5": []}}
    }


def test_recursive_dict_merge():
    defaults = {"a": 1, "b": {"a": 1, "b": 2}}
    calc_swaps = {"c": None, "b": {"b": 3, "d": 1}}
    assert recursive_dict_merge(defaults, calc_swaps) == {
        "a": 1,
        "b": {"a": 1, "b": 3, "d": 1},
        "c": None,
    }


def test_recursive_dict_merge2():
    defaults = {"a": 1, "b": {"a": 1, "b": 2}}
    calc_swaps = {"c": Remove, "b": {"b": 3, "d": 1}}
    assert recursive_dict_merge(defaults, calc_swaps) == {
        "a": 1,
        "b": {"a": 1, "b": 3, "d": 1},
    }


def test_recursive_dict_merge_verbose(caplog):
    defaults = {"a": 1, "b": {"a": 1, "b": 2}}
    calc_swaps = {"a": Remove, "b": {"b": 3, "d": 1}}
    with caplog.at_level(logging.WARNING):
        recursive_dict_merge(defaults, calc_swaps, verbose=True)
        assert "Overwriting key 'a'" in caplog.text
        assert "Overwriting key 'b' to: '3'" in caplog.text


def test_remove_instantiation():
    with pytest.raises(NotImplementedError):
        Remove()


def test_finalize_dict():
    with pytest.raises(ValueError, match="The directory should not"):
        finalize_dict({}, directory="tmp-quacc")
