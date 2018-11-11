BigDFT-suite manifesto
=======================

The code is not a monolithic piece of software but a collection of *independent*
packages that **may** be installed independently.
The :ref:`Installer <installer-script-label>` script has been
designed for the purpose.


Tips for developers
===================

Choose the correct package in which to insert the developments
--------------------------------------------------------------

The first question to ask yourself is the *generality* of the
functionality you are going to implement.
The spirit is to work at the lowest possible level for a given task.
The idea is to make available the functionality also to other potential
users of the ``BigDFT-suite`` subpackages.
This will also help in a better structure of the API of each package.

For instance, suppose you would like to implement a continuum solvent
cavity determination for a particular DFT run of a molecular system.
The correct level of development in this case would be the ``psolver``
package, as this is presently dealing with continuum solvents and cavities.

For a general overview one might say that:

 * ``futile`` deals with low-level functionalities like ``stdlib`` (but for FORTRAN).
   New MPI wrappers, strategies for memory copy and allocations should be implemented there.

 * ``at_lab`` library (will) deal with all the operations which are associated to position

 * <to_B_continued>


Read the coding rules
---------------------

For some inspiration on coding style and strategies, read this_.

.. _this: http://bigdft.org/Wiki/index.php?title=Coding_Rules


Document the API of the high-level routines
-------------------------------------------

.. todo:: write something here

Create a test for the functionality
-----------------------------------

Each of the packages has its own continuous integration procedure,  refer to
it for a suitable implementation.

 * ``futile``, ``psolver``: F_REGTEST_INSTRUCTION (to be documented)
 * ``bigdft`` see here_.

.. _here: http://bigdft.org/Wiki/index.php?title=Inserting_a_new_test_in_the_distribution

Make a notebook which demonstrates the functionality in PyBigDFT
----------------------------------------------------------------

For each new high level functionality, you should create a jupyter notebook which demonstrates the new capability.  The idea is to ensure continuity and to help acquaint users with the new feature.  Some examples of notebooks can be found on github_.

.. _github: https://github.com/luigigenovese/BigDFT-nb

Insert the notebook as a tutorial in the PyBigDFT documentation
---------------------------------------------------------------

Once an appropriate notebook has been written, this should be added to the tutorial directory (``BIGDFT_ROOT/PyBigDFT/source/tutorials``), so that the documentation will be automatically generated and available as a tutorial at :ref:`pybigdft:pybigdft_tutorials`.
