# Not To be installed:
noinst_LIBRARIES = liblibsiesta.a
module_DATA      = $(modules)


mod_src = m_region.F90 \
          intrinsic_missing.F90 \
    	  m_pivot_methods.F90 \
    	  class_Sparsity.F90 \
    	  class_OrbitalDistribution.F90 \
		  geom_helper.f90 \
		  m_die.F90

liblibsiesta_a_SOURCES = $(mod_src)

# Module files.
if CAPITALIZE
modules = M_REGION.@MODULE_EXT@ \
          INTRINSIC_MISSING.@MODULE_EXT@ \
    	  M_PIVOT_METHODS.@MODULE_EXT@ \
    	  CLASS_SPARSITY.@MODULE_EXT@ \
    	  CLASS_ORBITALDISTRIBUTION.@MODULE_EXT@ \
    	  GEOM_HELPER.@MODULE_EXT@ \
		  M_DIE.@MODULE_EXT@
else
modules = m_region.@MODULE_EXT@ \
          intrinsic_missing.@MODULE_EXT@ \
    	  m_pivot_methods.@MODULE_EXT@ \
    	  class_sparsity.@MODULE_EXT@ \
    	  class_orbitaldistribution.@MODULE_EXT@ \
    	  geom_helper.@MODULE_EXT@ \
		  m_die.@MODULE_EXT@

endif

CLEANFILES = *.@MODULE_EXT@

AM_FCFLAGS = -I. -I$(srcdir)  @LIB_FUTILE_CFLAGS@ @LIBETSFIO_INCLUDE@ @MPI_INCLUDE@

#dependencies
$(modules): $(liblibsiesta_a_OBJECTS)
m_pivot_methods.o: m_region.o
m_region.o: class_sparsity.o class_orbitaldistribution.o geom_helper.o
class_sparsity.o: m_die.o
class_orbitaldistribution.o: m_die.o
geom_helper.o: m_die.o
intrinsic_missing.o:
m_die.o:

# Copy modules in a common include directory.
all: all-am insmod
install: install-am insmod
check: check-am insmod
insmod: liblibsiesta.a
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done
