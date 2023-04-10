# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys

# import typing
# typing.TYPE_CHECKING = True
from quacc import __version__

sys.path.insert(0, os.path.abspath("../../"))

# -- Project information -----------------------------------------------------

project = "quacc"
copyright = "2021"
author = "Andrew S. Rosen"

# The short X.Y version
version = __version__
# The full version, including alpha/beta/rc tags
release = __version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "m2r2",
    "nbsphinx",
    "nbsphinx_link",
    "sphinxcontrib.autodoc_pydantic",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["Thumbs.db", ".DS_Store", "test*.py"]

# use type hints
autodoc_typehints = "description"
autoclass_content = "both"
autodoc_member_order = "bysource"

# better napolean support
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_use_ivar = True

# The suffix(es) of source filenames.
source_suffix = [".rst", ".md"]

# control pydantic model docs
autodoc_pydantic_model_show_json = False
autodoc_pydantic_model_show_field_summary = False
autodoc_pydantic_model_show_config = False
autodoc_pydantic_model_show_config_summary = False
autodoc_pydantic_model_show_validator_members = False
autodoc_pydantic_model_member_order = "bysource"
autodoc_pydantic_settings_show_json = False
autodoc_pydantic_settings_show_field_summary = False
autodoc_pydantic_settings_show_config = False
autodoc_pydantic_settings_show_config_summary = False
autodoc_pydantic_settings_show_validator_members = False
autodoc_pydantic_settings_member_order = "bysource"
autodoc_pydantic_field_list_validators = False
autodoc_pydantic_field_show_constraints = False

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "furo"

# hide sphinx footer
html_show_sphinx = False
html_show_sourcelink = False

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
fonts = [
    "Lato",
    "-apple-system",
    "BlinkMacSystemFont",
    "Segoe UI",
    "Helvetica",
    "Arial",
    "sans-serif",
    "Apple Color Emoji",
    "Segoe UI Emoji",
]
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_favicon = "_static/favicon.ico"
html_theme_options = {
    "light_css_variables": {
        "admonition-font-size": "92%",
        "admonition-title-font-size": "92%",
        "font-stack": ",".join(fonts),
        "font-size--small": "92%",
        "font-size--small--2": "87.5%",
        "font-size--small--3": "87.5%",
        "font-size--small--4": "87.5%",
    },
    "dark_css_variables": {
        "admonition-font-size": "92%",
        "admonition-title-font-size": "92%",
        "font-stack": ",".join(fonts),
        "font-size--small": "92%",
        "font-size--small--2": "87.5%",
        "font-size--small--3": "87.5%",
        "font-size--small--4": "87.5%",
    },
}
html_title = "quacc"

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    "python": ("https://docs.python.org/3.10", None),
    "ase": ("https://wiki.fysik.dtu.dk/ase/py-modindex.html", None),
    "pymatgen": ("http://pymatgen.org", None),
    "atomate2": ("https://materialsproject.github.io/atomate2/", None),
    "jobflow": ("https://materialsproject.github.io/jobflow/", None),
    "monty": ("https://guide.materialsvirtuallab.org/monty/", None),
    "cclib": ("https://cclib.github.io/", None),
    "emmet": ("https://materialsproject.github.io/emmet/", None),
}
