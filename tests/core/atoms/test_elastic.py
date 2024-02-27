from pathlib import Path

from ase.io import read

from quacc.atoms.deformation import make_deformations_from_bulk

FILE_DIR = Path(__file__).parent


def test_make_deformations_from_bulk():
    atoms = read(FILE_DIR / "ZnTe.cif.gz")
    atoms.info["test"] = "hi"
    deformations = make_deformations_from_bulk(atoms)
    assert len(deformations) == 24
    assert [
        deformation.get_atomic_numbers() == [30, 30, 30, 30, 52, 52, 52, 52]
        for deformation in deformations
    ]
    assert [
        deformation.get_chemical_formula() == "Te4Zn4" for deformation in deformations
    ]
    assert [atoms.info["test"] == "hi" for deformation in deformations]
