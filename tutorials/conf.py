# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import os
import sys

# Safely import project version
try:
    from melektrodica import __version__ as melektrodica_version

    version = melektrodica_version
    release = melektrodica_version
except ImportError:
    version = "unknown"
    release = "unknown"
    sys.stderr.write("Warning: 'μElektrodica' module not found. Unable to set version.\n")

# Add the root project directory to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../")))

# -- Project information -----------------------------------------------------

project = "μElektrodica"
copyright = "Copyright (C) 2024 C. Baqueiro Basto, M. Secanell, L.C. Ordoñez"
author = "μElektrodica Developers"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or custom ones.
language = 'en'

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "nbsphinx",
    "nbsphinx_link",
    "sphinxcontrib.bibtex",
]

# Sphinx-napoleon configuration (for Google/Numpy style docstrings)
napoleon_google_docstring = True

# Configuration for adding bibliographies
bibtex_bibfiles = ["bib.bib"]

# Paths for templates and ignored files
templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
    "**/temp_files/",
]

# -- Options for HTML output -------------------------------------------------

# Theme configuration
html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 3,
}

# Uncomment to customize styles
# html_static_path = ["_static"]
# html_style = "css/my_style.css"
