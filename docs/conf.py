# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

# -- Project information -----------------------------------------------------

project = 'EXCEED-DM'
copyright = '2022, Tanner Trickle'
author = 'Tanner Trickle'

# The full version, including alpha/beta/rc tags
version = '1.0.0'
release = '1.0.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
        'myst_parser',
        # 'sphinxfortran.fortran_domain',
        # 'sphinxfortran.fortran_autodoc',
        'sphinxcontrib.bibtex',
        'sphinx_copybutton',
        'sphinx_tabs.tabs'
]

# bibtex
bibtex_bibfiles = ['bibliography.bib']

bibtex_default_style = 'plain'

# Fortran
# fortran_src = [
#                     "../src/*.f90", 
#                     "../src/utils/*.f90"
#             ]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

html_theme = 'furo'
html_logo = "./media/exdm-prelim-logo-modified.png"
html_title = "EXCEED-DM v1.0.0"
html_favicon = "./media/exdm-favicon2.png"

# html_theme = 'sphinx_rtd_theme'

# html_theme = 'sphinx_book_theme'
# html_theme_options = {
#     "repository_url": "https://github.com/tanner-trickle/EXCEED-DM",
#     "use_repository_button": True,
#     "show_toc_level": 2,
#     "show_navbar_depth": 2,
# }

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

myst_enable_extensions = ["dollarmath", "smartquotes"]
