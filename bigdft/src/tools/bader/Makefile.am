bin_PROGRAMS = bader

bader_objects = kind_mod.f90       \
	matrix_mod.f90     \
	ions_mod.f90       \
	options_mod.f90    \
	charge_mod.f90     \
	chgcar_mod.f90     \
	cube_mod.f90       \
	io_mod.f90         \
	bader_mod.f90      \
	dipole_mod.f90     \
	voronoi_mod.f90

bader_SOURCES = main.f90 $(bader_objects)

#Use the flib library (use only yaml_output)
bader_LDADD = @LIB_FUTILE_LIBS@

AM_FCFLAGS = -I. @LIB_FUTILE_CFLAGS@
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS = -Wl,-rpath=$(DESTDIR)$(libdir)
endif

CLEANFILES = *.@MODULE_EXT@ 

#dependencies
main.o : $(bader_objects:.f90=.o)

ions_mod.o matrix_mod.o options_mod.o: kind_mod.o

charge_mod.o chgcar_mod.o bader_mod.o dipole_mod.o voronoi_mod.o : matrix_mod.o ions_mod.o

io_mod.o cube_mod.o : ions_mod.o charge_mod.o

dipole_mod.o chgcar_mod.o : charge_mod.o

bader_mod.o voronoi_mod.o dipole_mod.o: io_mod.o

io_mod.o : chgcar_mod.o cube_mod.o

dipole_mod.o : bader_mod.o
