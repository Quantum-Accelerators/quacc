from __future__ import annotations

from unittest.mock import Mock

import numpy as np
import pytest

from quacc.calculators.qchem.io import write_qchem


@pytest.mark.parametrize("coefficients", [[0.1, 0.2], np.array([0.1, 0.2])])
def test_write_orbital_coefficients(tmp_path, coefficients):
    qc_input = Mock()

    write_qchem(qc_input, tmp_path, prev_orbital_coeffs=coefficients)

    qc_input.write_file.assert_called_once_with(tmp_path / "mol.qin")
    np.testing.assert_array_equal(np.fromfile(tmp_path / "53.0"), coefficients)


@pytest.mark.parametrize("coefficients", [None, [], np.array([])])
def test_write_without_orbital_coefficients(tmp_path, coefficients):
    qc_input = Mock()

    write_qchem(qc_input, tmp_path, prev_orbital_coeffs=coefficients)

    qc_input.write_file.assert_called_once_with(tmp_path / "mol.qin")
    assert not (tmp_path / "53.0").exists()
