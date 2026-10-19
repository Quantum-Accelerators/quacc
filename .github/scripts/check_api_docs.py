"""Smoke-check the rendered public reference pages after building the docs."""

from __future__ import annotations

from html.parser import HTMLParser
from pathlib import Path


class ReferenceHeadings(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.ids: set[str] = set()

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag in {"h1", "h2", "h3", "h4", "h5", "h6"}:
            self.ids.update(value for key, value in attrs if key == "id" and value)


pages = list(Path("site/reference/quacc").rglob("*.html"))
assert pages, "Build the reference documentation first"
for page in pages:
    headings = ReferenceHeadings()
    headings.feed(page.read_text(encoding="utf-8"))
    assert not any(
        identifier.rsplit(".", 1)[-1].startswith("has_") for identifier in headings.ids
    ), page
for module in ("defects", "slabs"):
    page = Path(f"site/reference/quacc/recipes/common/{module}.html")
    html = page.read_text(encoding="utf-8")
    assert "ResultSchema" in html, page
    assert "quacc.types.RunSchema" in html, page
torchsim = Path("site/reference/quacc/recipes/torchsim/_base.html").read_text(
    encoding="utf-8"
)
for field in (
    "TrajectoryReporterDict.state_frequency",
    "AutobatcherDict.memory_scaling_factor",
):
    assert f"quacc.recipes.torchsim._base.{field}" in torchsim, field
print(f"Checked {len(pages)} rendered reference pages")  # noqa: T201
