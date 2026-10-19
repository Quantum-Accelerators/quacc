"""Generate compact schema diagrams without importing optional scientific packages."""

from __future__ import annotations

import ast
from html import escape
from pathlib import Path


def schema_diagram(source: str, roots: tuple[str, ...]) -> str:
    """Render schema fields, inheritance, and nested schemas from type annotations.

    Reading the syntax tree also handles classes inside TYPE_CHECKING blocks.
    External objects and open-ended dictionary aliases remain leaf field types.
    """
    classes = {
        node.name: node
        for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.ClassDef)
    }
    pending = list(roots)
    seen = set()
    nodes = []
    edges = []
    while pending:
        name = pending.pop(0)
        if name in seen:
            continue
        seen.add(name)
        cls = classes[name]
        optional = any(
            keyword.arg == "total" and ast.literal_eval(keyword.value) is False
            for keyword in cls.keywords
        )
        fields = []
        for field in cls.body:
            if not isinstance(field, ast.AnnAssign):
                continue
            field_name = ast.unparse(field.target)
            annotation = ast.unparse(field.annotation)
            marker = (
                " (optional key)"
                if optional or annotation.startswith("NotRequired[")
                else ""
            )
            fields.append(escape(f"{field_name}: {annotation}{marker}", quote=True))
            for child in sorted(
                {
                    node.id
                    for node in ast.walk(field.annotation)
                    if isinstance(node, ast.Name)
                    and node.id.endswith("Schema")
                    and node.id in classes
                }
            ):
                pending.append(child)
                edges.append(f'    {name} -->|"{field_name}"| {child}')
        label = "<br/>".join([f"<b>{name}</b>", *fields])
        nodes.append(f'    {name}["{label}"]')
        for base in cls.bases:
            if isinstance(base, ast.Name) and base.id in classes:
                pending.append(base.id)
                edges.append(f'    {name} -.->|"inherits fields"| {base.id}')

    links = " · ".join(f"[{name}][quacc.types.{name}]" for name in sorted(seen))
    return (
        "```mermaid\nflowchart TB\n"
        + "\n".join(nodes + edges)
        + "\n```\n\nDotted arrows mean inherited fields; solid arrows mean nested values. "
        "External objects and calculator-specific dictionaries are shown as field types.\n\n"
        + links
        + "\n"
    )


def write_schema_diagrams(out_dir: Path) -> None:
    """Write snippets consumed by recipe docstrings during the docs build."""
    source = Path("src/quacc/types.py").read_text(encoding="utf-8")
    target = out_dir / "_schemas"
    target.mkdir(parents=True, exist_ok=True)
    (target / "ase.md").write_text(
        schema_diagram(source, ("OptSchema", "RunSchema")), encoding="utf-8"
    )
