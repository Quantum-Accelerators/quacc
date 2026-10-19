# Updating the Docs

Install the documentation dependencies with `pip install -e ".[docs]"`.

The documentation uses [Zensical](https://zensical.org/), configured in `_zensical.toml`. From the repository root, generate the reference pages and navigation before building:

```bash
python .github/scripts/gen_ref_pages.py
python .github/scripts/gen_zensical_toml.py
zensical build --clean
python .github/scripts/check_api_docs.py
python -m http.server --directory site 8000
```

Open `http://localhost:8000` to preview the result. Check the affected reference pages for readable signatures, useful parameter and return descriptions, and working schema links. The reference renderer reads the local `src` tree and hides internal dependency flags and type variables. Generated reference pages, `zensical.toml`, and `site/` should not be committed.
