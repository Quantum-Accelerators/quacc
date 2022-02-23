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
copyright = "2022"
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
    "numpydoc",
    # "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "m2r2",
    "sphinx_panels",
    "sphinxcontrib.autodoc_pydantic",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["Thumbs.db", ".DS_Store", "test*.py"]

# use type hints
# autodoc_typehints = "description"
autoclass_content = "class"
autodoc_member_order = "bysource"
# autodoc_inherit_docstrings = True
python_use_unqualified_type_names = True

# autodoc_typehints_description_target = "documented"
# autodoc_class_signature = "separated"

autodoc_preserve_defaults = True

# don't overwrite summary but generate them if they don't exist
autosummary_generate = True
autosummary_generate_overwrite = True

# numpydoc options
numpydoc_class_members_toctree = False
numpydoc_show_class_members = False
numpydoc_show_inherited_class_members = False
numpydoc_attributes_as_param_list = False
numpydoc_xref_param_type = True

# better napolean support
# napoleon_use_param = True
# napoleon_use_rtype = True
# napoleon_use_ivar = True

# sphinx-panels shouldn't add bootstrap css as the pydata-sphinx-theme already loads it
panels_add_bootstrap_css = False

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

autosummary_imported_members = False

# The suffix(es) of source filenames.
source_suffix = [".rst", ".md"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "pydata_sphinx_theme"

# hide sphinx footer
html_show_sphinx = False
html_show_sourcelink = False
html_theme_options = {
    "github_url": "https://github.com/arosen93/quacc",
    # "use_edit_page_button": True,
    "show_toc_level": 1,
}

html_context = {
    "github_user": "arosen93",
    "github_repo": "quacc",
    "github_version": "main",
    "doc_path": "docs",
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = ["custom.css"]
# html_favicon = "_static/favicon.ico"

html_title = "quacc"

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    "python": ("https://docs.python.org/3.10", None),
    "matplotlib": ("http://matplotlib.org", None),
    "ase": ("https://wiki.fysik.dtu.dk/ase/index.html", None),
    "pymatgen": ("http://pymatgen.org", None),
    "atomate2": ("https://materialsproject.github.io/atomate2/", None),
    "jobflow": ("https://materialsproject.github.io/jobflow/", None),
    "monty": ("https://guide.materialsvirtuallab.org/monty/", None),
    "cclib": ("https://cclib.github.io/", None),
}
