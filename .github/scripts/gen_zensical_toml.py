"""
Using _zensical.toml as a starting file, generate a zensical.toml with navigation links
generated from `docs/reference/SUMMARY.md`.
This is roughly equivalent to what the `mkdocs-literate-nav` plugin would have done.
"""

from __future__ import annotations

import re
from pathlib import Path


def parse_summary(path: Path, prefix: str = "reference/") -> list:
    nav = []
    stack = [(nav, -1)]
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip())
        m = re.match(r"\s*\*\s+\[(.+?)\]\((.+?)\)", line)
        if not m:
            continue
        title, file_path = m.groups()

        while stack[-1][1] >= indent:
            stack.pop()

        parent_list, _ = stack[-1]
        children = []
        parent_list.append((title, prefix + file_path, children))
        stack.append((children, indent))
    return nav


def to_toml_nav(entries: list, indent: int = 8) -> str:
    pad = " " * indent
    lines = []
    for title, path, children in entries:
        if children:
            inner = to_toml_nav(children, indent + 4)
            lines.append(f'{pad}{{ "{title}" = [\n{inner}{pad}] }}')
        else:
            lines.append(f'{pad}{{ "{title}" = "{path}" }}')
    return ",\n".join(lines) + "\n"


summary = Path("docs/reference/SUMMARY.md")
nav_entries = parse_summary(summary)
inner_toml = to_toml_nav(nav_entries, indent=12)
replacement = f'{{ "Code Documentation" = [\n        {{ "quacc" = [\n{inner_toml}        ] }}\n    ] }}'

toml_path = Path("_zensical.toml")
content = toml_path.read_text()
content = re.sub(
    r'\{\s*"Code Documentation"\s*=\s*"reference/"\s*\}', replacement, content
)
Path("zensical.toml").write_text(content)
