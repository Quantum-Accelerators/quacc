""" "
Using _zensical.toml as a starting file, generate a zensical.toml with navigation links
generated from `docs/reference/SUMMARY.md`.

This is roughly equivalent to what the `mkdocs-literate-nav` plugin would have done.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


def parse_summary(path: Path) -> list[tuple[int, str, str | None]]:
    """Parse SUMMARY.md into a flat list of (indent_level, label, path_or_None)."""
    items = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        n_spaces = len(line) - len(line.lstrip(" "))
        indent = n_spaces // 4
        rest = line.strip()
        if not rest.startswith("* "):
            raise ValueError(f"Unexpected line format: {line!r}")
        rest = rest[2:]
        m = re.match(r"\[([^\]]+)\]\(([^)]+)\)", rest)
        if m:
            items.append((indent, m.group(1), m.group(2)))
        else:
            items.append((indent, rest, None))
    return items


def build_tree(items: list[tuple[int, str, str | None]]) -> list[dict]:
    """Build a nested tree from flat (indent, label, path) items."""
    root: list[dict] = []
    # stack entries: (indent_level, children_list)
    stack: list[tuple[int, list]] = [(-1, root)]

    for indent, label, path in items:
        node: dict = {"label": label, "path": path, "children": []}
        while stack[-1][0] >= indent:
            stack.pop()
        stack[-1][1].append(node)
        stack.append((indent, node["children"]))

    return root


def render_children(children: list[dict], indent: int, path_prefix: str) -> list[str]:
    """Render a list of nav nodes as indented TOML lines."""
    lines = []
    pad = "    " * indent
    for i, child in enumerate(children):
        comma = "," if i < len(children) - 1 else ""
        if not child["children"]:
            path = path_prefix + child["path"]
            lines.append(f'{pad}{{ "{child["label"]}" = "{path}" }}{comma}')
        else:
            lines.append(f'{pad}{{ "{child["label"]}" = [')
            lines.extend(render_children(child["children"], indent + 1, path_prefix))
            lines.append(f"{pad}] }}{comma}")
    return lines


def generate_code_docs_entry(
    tree: list[dict], indent: int = 1, path_prefix: str = "reference/"
) -> str:
    """Generate the full 'Code Documentation' nav entry as a TOML string."""
    pad = "    " * indent
    lines = [f'{pad}{{ "Code Documentation" = [']
    lines.extend(render_children(tree, indent + 1, path_prefix))
    lines.append(f"{pad}] }},")
    return "\n".join(lines)


def main() -> None:
    args = sys.argv[1:]
    summary_path = Path(args[0] if args else "docs/reference/SUMMARY.md")
    input_toml = Path(args[1] if len(args) > 1 else "_zensical.toml")
    output_toml = Path(args[2] if len(args) > 2 else "zensical.toml")

    items = parse_summary(summary_path)
    tree = build_tree(items)
    code_docs_entry = generate_code_docs_entry(tree)

    toml_content = input_toml.read_text()

    placeholder_re = re.compile(
        r'^[ \t]*\{ "Code Documentation" = "[^"]*" \},?[ \t]*$', re.MULTILINE
    )
    match = placeholder_re.search(toml_content)
    if not match:
        sys.exit(1)

    result = (
        toml_content[: match.start()]
        + code_docs_entry
        + "\n"
        + toml_content[match.end() + 1 :]
    )
    output_toml.write_text(result)


if __name__ == "__main__":
    main()
