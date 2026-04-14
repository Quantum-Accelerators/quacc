"""
Generate docs/reference/ files to disk.
This is roughly equivalent to what the `mkdocs-gen-files` plugin would have done.
"""

from __future__ import annotations

from pathlib import Path

out_dir = Path("docs/reference")
out_dir.mkdir(parents=True, exist_ok=True)

nav_lines = []
seen_dirs: set[tuple[str, ...]] = set()

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

    # Emit section headers for intermediate directories not yet seen
    for i in range(1, len(parts)):
        dir_key = parts[:i]
        if dir_key not in seen_dirs:
            seen_dirs.add(dir_key)
            indent = "    " * (len(dir_key) - 1)
            nav_lines.append(f"{indent}* {dir_key[-1]}\n")

    # Emit leaf link
    name = parts[-1]
    indent = "    " * (len(parts) - 1)
    nav_lines.append(f"{indent}* [{name}]({doc_path.as_posix()})\n")

(out_dir / "SUMMARY.md").write_text("".join(nav_lines))
