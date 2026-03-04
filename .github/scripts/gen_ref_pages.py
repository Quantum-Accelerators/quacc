"""
Generate docs/reference/ files to disk.
This is roughly equivalent to what the `mkdocs-gen-files` plugin would have done.
"""

from __future__ import annotations

from pathlib import Path

out_dir = Path("docs/reference")
out_dir.mkdir(parents=True, exist_ok=True)

nav_lines = []

for path in sorted(Path("src").rglob("*.py")):
    module_path = path.relative_to("src").with_suffix("")
    doc_path = path.relative_to("src").with_suffix(".md")
    full_doc_path = out_dir / doc_path

    parts = tuple(module_path.parts)
    ignore = ["_cli", "_version", "__init__", "__main__"]
    if any(p in ignore for p in parts):
        continue

    full_doc_path.parent.mkdir(parents=True, exist_ok=True)
    full_doc_path.write_text(f"::: {'.'.join(parts)}\n")

    indent = "    " * (len(parts) - 1)
    nav_lines.append(f"{indent}* [{parts[-1]}]({doc_path.as_posix()})\n")

(out_dir / "SUMMARY.md").write_text("".join(nav_lines))
