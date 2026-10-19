"""Tests for the dependency-free documentation schema generator."""

from __future__ import annotations

import unittest
from pathlib import Path

from schema_diagrams import schema_diagram


class SchemaDiagramTests(unittest.TestCase):
    def test_type_checking_inheritance_and_nested_containers(self):
        source = """
if TYPE_CHECKING:
    class OptionalSchema(TypedDict, total=False):
        metadata: dict[str, Any]
    class BaseSchema(OptionalSchema):
        atoms: Atoms
    class ResultSchema(BaseSchema):
        children: list[BaseSchema | None]
        note: NotRequired[str]
"""
        diagram = schema_diagram(source, ("ResultSchema",))
        assert 'ResultSchema -.->|"inherits fields"| BaseSchema' in diagram
        assert 'ResultSchema -->|"children"| BaseSchema' in diagram
        assert "metadata: dict[str, Any] (optional key)" in diagram
        assert "note: NotRequired[str] (optional key)" in diagram
        assert diagram.count('BaseSchema["') == 1
        assert 'Atoms["' not in diagram

    def test_cycles_and_html_escaping(self):
        source = """
class NodeSchema(TypedDict):
    children: list[NodeSchema]
    tag: Literal["<tag>"]
"""
        diagram = schema_diagram(source, ("NodeSchema",))
        assert diagram.count('NodeSchema["') == 1
        assert "&lt;tag&gt;" in diagram

    def test_current_ase_schemas(self):
        diagram = schema_diagram(
            Path("src/quacc/types.py").read_text(encoding="utf-8"),
            ("OptSchema", "RunSchema"),
        )
        assert "trajectory_results: list[Results]" in diagram
        assert "input_atoms: AtomsSchema | None" in diagram
        assert "structure_metadata: StructureMetadata (optional key)" in diagram
        assert "[OptSchema][quacc.types.OptSchema]" in diagram


if __name__ == "__main__":
    unittest.main()
