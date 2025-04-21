from __future__ import annotations

import pytest

from quacc.utils.kpts import bandgap_to_kspacing


def test_bandgap_to_kspacing():
    out = bandgap_to_kspacing(0.60)
    assert pytest.approx(0.2004)
