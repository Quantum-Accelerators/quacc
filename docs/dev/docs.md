# Updating the Docs

To install the dependencies to build the documentation, run `pip install -e .[docs]`.

The quacc documentation is built using [Zensical](https://zensical.org/), configured via `_zensical.toml` in the base directory. To build the documentation locally (e.g. to evaluate changes), run `zensical build --clean` in the base directory and open `site/index.html` in a browser.
