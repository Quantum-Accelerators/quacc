"""Regression tests for generated reference-page edit links."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path


def test_reference_edit_links(tmp_path):
    repo = Path(__file__).resolve().parents[2]
    source = tmp_path / "src/quacc/recipes/vasp/mp_legacy.py"
    source.parent.mkdir(parents=True)
    source.write_text('"""Example module."""\n')
    source.with_name("__init__.py").write_text("")
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "index.md").write_text("# Home\n")
    overrides = (repo / "docs/overrides").as_posix()
    config = (
        '[project]\nsite_name = "Edit link regression"\n'
        'repo_url = "https://github.com/Quantum-Accelerators/quacc/"\n'
        'edit_uri = "blob/main/docs/"\nuse_directory_urls = false\n'
        f'[project.theme]\ncustom_dir = "{overrides}"\n'
        'features = ["content.action.edit"]\n'
    )
    (tmp_path / "_zensical.toml").write_text(config)
    (tmp_path / "zensical.toml").write_text(config)
    subprocess.run(
        [sys.executable, str(repo / ".github/scripts/gen_ref_pages.py")],
        cwd=tmp_path,
        check=True,
    )
    reference = docs / "reference/quacc/recipes/vasp/mp_legacy.md"
    expected = (
        "https://github.com/Quantum-Accelerators/quacc/"
        "blob/main/src/quacc/recipes/vasp/mp_legacy.py"
    )
    assert f"edit_url: {expected}" in reference.read_text()
    assert "::: quacc.recipes.vasp.mp_legacy" in reference.read_text()
    assert not reference.with_name("__init__.md").exists()
    subprocess.run(
        [sys.executable, "-m", "zensical", "build", "--clean"], cwd=tmp_path, check=True
    )
    for page, url in (
        ("reference/quacc/recipes/vasp/mp_legacy.html", expected),
        (
            "index.html",
            "https://github.com/Quantum-Accelerators/quacc//blob/main/docs/index.md",
        ),
    ):
        html = (tmp_path / "site" / page).read_text(encoding="utf-8")
        edit_links = re.findall(r'<a\b[^>]*href="([^"]+)"[^>]*rel="edit"', html)
        assert edit_links == [url]
