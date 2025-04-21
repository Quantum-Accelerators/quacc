from __future__ import annotations

import math

from quacc.utils.kpts import bandgap_to_kspacing


def test_bandgap_to_kspacing():
    out = bandgap_to_kspacing(0.60)
    assert math.isclose(out, 0.2004, rel_tol=1e-6)
