# To be installed:
dynlibdir        = $(libdir)
lib_LIBRARIES    = libatlab-1.a
dynlib_DATA      = $(atlab_dynlib) $(openbabel_dynamic_library)
noinst_DATA      = $(atlab_linklib) $(openbabel_soname_library)
module_DATA      = $(modules)

# Static or dynamic library
atlab_dynamic_library = libatlab-1.so.@ATLAB_MINOR_VERSION@.0.@ATLAB_MICRO_VERSION@
atlab_library_soname  = libatlab-1.so.@ATLAB_MINOR_VERSION@

# We install atlab library only if PSolver is deactivated.
if BUILD_DYNAMIC_LIBS
atlab_dynlib = $(atlab_dynamic_library)
atlab_linklib = $(atlab_library_soname)
else
atlab_dynlib =
atlab_linklib =
endif

$(modules): $(libatlab_1_a_OBJECTS)

libatlab_1_a_SOURCES = f_functions.f90 \
    harmonics.f90 \
    box.f90 \
    domain.f90  \
    IObox.f90 \
    ISF.f90 \
    mp_quadrature.f90 \
    numerics.f90 \
    $(openbabel_static_sources) \
    openbabel_wrapper.f90 \
    $(etsf_files)


EXTRA_DIST= openbabel.cpp FiniteDiffCorff.inc

f_functions.o: numerics.o box.o FiniteDiffCorff.inc
harmonics.o: numerics.o box.o
box.o: domain.o
domain.o: numerics.o
mp_quadrature.o IObox.o IOetsf.o IOetsf_fake.o: box.o

if HAVE_ETSF_IO
etsf_files = IOetsf.f90
etsf_dep   = IOetsf.o
else
etsf_files = IOetsf_fake.f90
etsf_dep   = IOetsf_fake.o
endif

AM_FCFLAGS = -I. -I$(srcdir) @LIB_FUTILE_CFLAGS@ @LIBETSFIO_INCLUDE@

if HAVE_OPENBABEL
openbabel_dynamic_sources = openbabel.cpp
openbabel_dynamic_library = libbabel.so.@ATLAB_MINOR_VERSION@.0.@ATLAB_MICRO_VERSION@
openbabel_soname_library = libbabel.so.@ATLAB_MINOR_VERSION@
openbabel_library = $(openbabel_soname_library)
AM_LDFLAGS = @FC_RPATH@$(DESTDIR)$(dynlibdir)
else
openbabel_static_sources = openbabel_fake.c
AM_CFLAGS = -I. -I$(srcdir) -I$(top_builddir)/includes @LIB_FUTILE_CFLAGS@
endif

AM_CXXFLAGS = @LIB_FUTILE_CFLAGS@ @LIB_OPENBABEL_CFLAGS@

# Module files.
if CAPITALIZE
modules = F_FUNCTIONS.@MODULE_EXT@ \
	F_HARMONICS.@MODULE_EXT@ \
	NUMERICS.@MODULE_EXT@ \
	AT_DOMAIN.@MODULE_EXT@ \
	AT_BABEL.@MODULE_EXT@ \
	BOX.@MODULE_EXT@ \
	IOBOXETSF.@MODULE_EXT@ \
	MULTIPOLE_PRESERVING.@MODULE_EXT@ \
	IOBOX.@MODULE_EXT@
else
modules = f_functions.@MODULE_EXT@ \
	f_harmonics.@MODULE_EXT@ \
	numerics.@MODULE_EXT@ \
	at_domain.@MODULE_EXT@ \
	at_babel.@MODULE_EXT@ \
	box.@MODULE_EXT@ \
	ioboxetsf.@MODULE_EXT@ \
	multipole_preserving.@MODULE_EXT@ \
	iobox.@MODULE_EXT@
endif
CLEANFILES = *.@MODULE_EXT@ libatlab-1.so* libbabel.so*

IObox.o: $(etsf_dep) harmonics.o

# Dynamic library building
$(atlab_dynamic_library): $(libatlab_1_a_OBJECTS) $(openbabel_library)
	$(FC) $(FCFLAGS) @FC_BUILD_SHARED@ $(LDFLAGS) @FC_RPATH@$(DESTDIR)$(dynlibdir) -Wl,-E -Wl,-soname=$(atlab_library_soname) -o $@ $(libatlab_1_a_OBJECTS) @LIBS_EXT_DEPS@ $(LIBS) $(openbabel_library)
	@chmod a+x $@

$(atlab_library_soname): $(atlab_dynamic_library)
	ln -fs $^ $@
	ln -fs $@ libatlab-1.so

$(openbabel_dynamic_library): $(openbabel_dynamic_sources)
	$(CXX) $(CXXFLAGS) @CC_PIC_FLAG@ @LIB_OPENBABEL_CFLAGS@ -shared -I$(top_builddir) @GLIB_CFLAGS@ @LIB_FUTILE_CFLAGS@ $(LDFLAGS) $(AM_LDFLAGS) -Wl,-soname=$(openbabel_soname_library) -o $@ $< @LIB_OPENBABEL_LIBS@ @LIB_LIBDICTS_LIBS@
	@chmod a+x $@

$(openbabel_soname_library): $(openbabel_dynamic_library)
	ln -fs $^ $@
	ln -fs $@ libbabel.so

install-data-hook:
	if test -n "$(openbabel_dynamic_library)" ; then \
	  cd $(DESTDIR)$(dynlibdir) ; \
	  chmod a+x $(openbabel_dynamic_library) ; \
	  ln -fs $(openbabel_dynamic_library) $(openbabel_soname_library) ; \
	  ln -fs $(openbabel_soname_library) libbabel.so ; \
	fi
	if test -n "$(atlab_dynlib)" ; then \
	  cd $(DESTDIR)$(dynlibdir) ; \
	  chmod a+x $(atlab_dynamic_library) ; \
	  ln -fs $(atlab_dynamic_library) $(atlab_library_soname) ; \
	  ln -fs $(atlab_library_soname) libatlab-1.so ; \
	fi
