# Updating the Docs

To install the dependencies to build the documentation, run `pip install -e .[docs]`.

The quacc documentation is built using [Zensical](https://zensical.org/), configured via `_zensical.toml` in the base directory. To build the documentation locally (e.g. to evaluate changes), run `zensical build --clean` in the base directory and open `site/index.html` in a browser.

First generate the reference pages (including schema diagrams) and configuration,
as in the documentation CI workflow:

```bash
python .github/scripts/gen_ref_pages.py
python .github/scripts/gen_zensical_toml.py
zensical build --clean
```

The EMT `bulk_to_slabs_flow` return documentation includes an expandable output
diagram. `.github/scripts/schema_diagrams.py` generates its fields and relationships
directly from `src/quacc/types.py`, including classes defined under `TYPE_CHECKING`.
It uses only the Python standard library; no calculator imports or Graphviz
installation are needed. Dotted arrows show inheritance, solid arrows show nested
schemas, and optional keys are labeled. External objects and arbitrary calculator
dictionaries stay compact. The recipe docstring describes the flow-specific list
structure and includes the generated Markdown snippet.
