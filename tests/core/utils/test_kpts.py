from __future__ import annotations

import pytest

from quacc.utils.kpts import bandgap_to_kspacing


def test_bandgap_to_kspacing():
    assert bandgap_to_kspacing(0.60) == pytest.approx(0.2004)
