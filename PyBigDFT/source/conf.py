# -*- coding: utf-8 -*-
#
# PyBigDFT documentation build configuration file, created by
# sphinx-quickstart on Thu Apr 19 10:13:37 2018.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
from os import path as p
import sys
sourcepath=p.abspath(p.realpath(__file__))
pybigdftpath=p.dirname(p.dirname(sourcepath))
bigdftsuitepath=p.dirname(pybigdftpath)
sys.path.insert(0, p.abspath(p.join(bigdftsuitepath,'futile','src','python')))
sys.path.insert(0, p.abspath(pybigdftpath))
sys.path.insert(0, p.abspath(p.join(pybigdftpath,'notebooks')))

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    #'sphinxcontrib.jinja',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary']#,
  #'nbsphinx']

# import guzzle_sphinx_theme

# html_theme_path = guzzle_sphinx_theme.html_theme_path()
# html_theme = 'guzzle_sphinx_theme'

# Register the theme as an extension to generate a sitemap.xml
# extensions.append("guzzle_sphinx_theme")

# Guzzle theme options (see theme.conf for more information)
#html_theme_options = {
#    # Set the name of the project to appear in the sidebar
#    "project_nav_name": "PyBigDFT",
#}


# Add any paths that contain templates here, relative to this directory.
templates_path = ['sphinx_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'PyBigDFT'
copyright = u'2018, Luigi Genovese'
author = u'Luigi Genovese'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = u'0.1'
# The full version, including alpha/beta/rc tags.
release = u'0.1'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'bizstyle'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['sphinx_static']
html_logo = 'sphinx_static/logo_header.png'
html_favicon = 'sphinx_static/logo_header.png'


# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'PyBigDFTdoc'


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'PyBigDFT.tex', u'PyBigDFT Documentation',
     u'Luigi Genovese', 'manual'),
]


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'pybigdft', u'PyBigDFT Documentation',
     [author], 1)
]


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'PyBigDFT', u'PyBigDFT Documentation',
     author, 'PyBigDFT', 'One line description of project.',
     'Miscellaneous'),
]



# -- Options for Epub output ----------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project
epub_author = author
epub_publisher = author
epub_copyright = copyright

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']



# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'https://docs.python.org/': None}

from BigDFT import InputActions as A

# information about the input actions that can be used and their documentation
# jinja_contexts = {'input_actions': { 'actions': { a: getattr(A,a).__doc__  for a in dir(A) if '__' not in a and a !='dict_set'}}}

# tools to be done in the readthedocs environment
#import os
#on_rtd = os.environ.get('READTHEDOCS') == 'True'
#if on_rtd:
#    sys.path.insert(0, p.abspath('../../futile/src/python/'))
#    import sys
#    from unittest.mock import MagicMock
#
#    class Mock(MagicMock):
#        @classmethod
#        def __getattr__(cls, name):
#            return MagicMock()
#
#    MOCK_MODULES = ['yaml','gi.repository']
#    sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)
