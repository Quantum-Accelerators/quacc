from atomate2.common.schemas.structure import StructureMetadata
from pymatgen.io.ase import AseAtomsAdaptor


def atoms_to_db(atoms, **kwargs):
    struct = AseAtomsAdaptor().get_structure(atoms)
    struct_metadata = StructureMetadata().from_structure(struct, **kwargs)
    return struct_metadata
