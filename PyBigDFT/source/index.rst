.. PyBigDFT documentation master file, created by
   sphinx-quickstart on Thu Apr 19 10:13:37 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PyBigDFT's documentation!
====================================

This package is a collection of Python Modules that are conceived for pre- and post- processing
of BigDFT input and output files. Such modules are supposed to enhance the BigDFT experience
by high-level approach.
Also, calculators and workflows are supposed to be created and inspected with modules of the PyBigDFT package.

How to write documentation in this package
==========================================

The syntax in the text pages is ``ReStructured Text`` format.
see `this page`_ for examples of the syntax, and `this one`__ for a more comprehensive explanation of
the syntax.

.. __: http://docutils.sourceforge.net/rst.html

.. _this page: http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html

The syntax chosen for the auto-documentation of the Python module is ``Google Docstring`` syntax.
see `the google Docstring example`__ in the `Sphinx-napoleon`__ package.

.. __: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html

.. __: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/index.html

How to generate documentation
-----------------------------

 * Install ``sphinx`` package and the requirements written in the file  ``source/requirements.txt`` of
   provided with the package.

 * Copy the quick build file ``Makefile-sphinxbuild`` of the ``PyBigDFT`` source directory into ``Makefile``

 * Run the command ``make html``

 * The documentation of the package can be explored starting from the file ``build/html/index.html``.

Automodule generation
---------------------

The directive to be used is mainly the ``automodule`` one, see `this`__ page:

.. __: http://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html

Module Members
==============

Input file representation  and usage from Python

.. toctree:: inputfiles

Input Actions on the BigDFT input file

.. toctree:: inputactions

Analysis of Logfiles and Ground-state related properties

.. toctree:: BigDFT.Logfiles

Processing of Fragment-related quantities

.. toctree:: BigDFT.Fragments

Calculators to be used with BigDFT package to trigger runs of the code

.. toctree:: BigDFT.Calculators

Organize runs and analyze output in a dataset

.. toctree:: datasets

(Py)BigDFT Tutorial page

.. toctree::
   :maxdepth: 1

   tutorials

Index of the Package
====================

.. toctree:: modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
