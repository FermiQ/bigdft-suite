# Not To be installed:
noinst_LIBRARIES = libsparsematrix.a
module_DATA      = $(modules)

# Code source.
if USE_BLACS
BLACS =
else
BLACS = blacs_fake.f90
endif

mod_src = sparsematrix_base.f90 \
          sparsematrix_init.f90 \
          sparsematrix.f90 \
    	  sparsematrix_io.f90 \
    	  parallel_linalg.f90 \
    	  sparsematrix_errorhandling.f90 \
    	  sparsematrix_timing.f90 \
    	  sparsematrix_types.f90 \
    	  sparsematrix_memory.f90 \
    	  sparsematrix_defs.f90 \
    	  m_region.F90 \
          intrinsic_missing.F90 \
    	  m_pivot_methods.F90 \
    	  class_Sparsity.F90 \
    	  class_OrbitalDistribution.F90 \
		  geom_helper.f90 \
		  m_die.F90 \
    	  $(BLACS)

libsparsematrix_a_SOURCES = $(mod_src)

# Module files.
if CAPITALIZE
modules = SPARSEMATRIX_BASE.@MODULE_EXT@ \
          SPARSEMATRIX_INIT.@MODULE_EXT@ \
          SPARSEMATRIX.@MODULE_EXT@ \
    	  PARALLEL_LINALG.@MODULE_EXT@ \
    	  SPARSEMATRIX_IO.@MODULE_EXT@ \
    	  SPARSEMATRIX_ERRORHANDLING.@MODULE_EXT@ \
    	  SPARSEMATRIX_TIMING.@MODULE_EXT@ \
    	  SPARSEMATRIX_TYPES.@MODULE_EXT@ \
    	  SPARSEMATRIX_MEMORY.@MODULE_EXT@ \
    	  SPARSEMATRIX_DEFS.@MODULE_EXT@ \
    	  M_REGION.@MODULE_EXT@ \
          INTRINSIC_MISSING.@MODULE_EXT@ \
    	  M_PIVOT_METHODS.@MODULE_EXT@ \
    	  CLASS_SPARSITY.@MODULE_EXT@ \
    	  CLASS_ORBITALDISTRIBUTION.@MODULE_EXT@ \
    	  GEOM_HELPER.@MODULE_EXT@ \
		  M_DIE.@MODULE_EXT@
else
modules = sparsematrix_base.@MODULE_EXT@ \
          sparsematrix_init.@MODULE_EXT@ \
          sparsematrix.@MODULE_EXT@ \
    	  parallel_linalg.@MODULE_EXT@ \
    	  sparsematrix_io.@MODULE_EXT@ \
    	  sparsematrix_errorhandling.@MODULE_EXT@ \
    	  sparsematrix_timing.@MODULE_EXT@ \
    	  sparsematrix_types.@MODULE_EXT@ \
    	  sparsematrix_memory.@MODULE_EXT@ \
    	  sparsematrix_defs.@MODULE_EXT@ \
    	  m_region.@MODULE_EXT@ \
          intrinsic_missing.@MODULE_EXT@ \
    	  m_pivot_methods.@MODULE_EXT@ \
    	  class_sparsity.@MODULE_EXT@ \
    	  class_orbitaldistribution.@MODULE_EXT@ \
    	  geom_helper.@MODULE_EXT@ \
		  m_die.@MODULE_EXT@
endif

CLEANFILES = *.@MODULE_EXT@

AM_FCFLAGS = -I. -I$(srcdir) @LIB_FUTILE_CFLAGS@ @LIBETSFIO_INCLUDE@ @MPI_INCLUDE@

#dependencies
$(modules): $(libsparsematrix_a_OBJECTS)
sparsematrix_defs.o: 
sparsematrix_memory.o: sparsematrix_errorhandling.o sparsematrix_types.o sparsematrix_defs.o
sparsematrix_base.o: sparsematrix_errorhandling.o sparsematrix_timing.o sparsematrix_types.o sparsematrix_memory.o sparsematrix_defs.o
parallel_linalg.o: sparsematrix_base.o sparsematrix_init.o
sparsematrix_init.o: sparsematrix_base.o m_pivot_methods.o
sparsematrix.o: sparsematrix_base.o sparsematrix_init.o parallel_linalg.o
sparsematrix_io.o: sparsematrix_base.o sparsematrix.o
m_pivot_methods.o: m_region.o
m_region.o: intrinsic_missing.o class_sparsity.o class_orbitaldistribution.o geom_helper.o
intrinsic_missing.o:
class_sparsity.o: m_die.o
class_orbitaldistribution.o: m_die.o
geom_helper.o: m_die.o
m_die.o:

# Copy modules in a common include directory.
all: all-am insmod
install: install-am insmod
check: check-am insmod
insmod: libsparsematrix.a
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done
