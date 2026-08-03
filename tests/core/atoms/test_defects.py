from __future__ import annotations

import pytest

pytest.importorskip("pymatgen.analysis.defects")

from pymatgen.core import Lattice, Structure

from quacc.atoms.defects import get_defect_entry_from_defect


def test_get_defect_entry_requires_dummy_site():
    defect_supercell = Structure(Lattice.cubic(5), ["Si"], [[0, 0, 0]])

    with pytest.raises(ValueError, match="No dummy site found in defect supercell"):
        get_defect_entry_from_defect(object(), defect_supercell, 0)
