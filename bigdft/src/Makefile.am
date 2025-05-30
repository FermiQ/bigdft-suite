
## Process this file with automake to produce Makefile.in

# Main targets.
# ------------
dynlibdir           = $(libdir)
bin_PROGRAMS        = $(binaries_sources)
noinst_PROGRAMS     = $(local_binaries)
bin_SCRIPTS         = $(scripts_sources)
lib_LIBRARIES       = libbigdft-1.a
module_DATA         = $(mod_bigdft)
dynlib_DATA         = $(bigdft_dynamic_library)

# Sub-directories.
# ---------------
SUBDIRS = \
	locregs \
	orbitals \
	modules \
	xcPSolver \
	$(sub_dirs_libconfig) \
	CUDA \
	OpenCL \
	bindings \
	geopt \
	$(sub_dirs_bin)

EXTRA_DIST = \
	wavelib/intots.inc \
	wavelib/recs16.inc \
	wavelib/sym_16.inc \
	init/pspconf.in.f90

# Build targets, binaries and libraries.
# -------------
noinst_LIBRARIES = libmain.a libabscalc.a
# put all subs
subs_static_library = \
	modules/libmodules.a \
	CUDA/libCUDA.a \
	locregs/liblocregs.a \
	orbitals/liborbitals.a \
	geopt/libgeopt.a \
	xcPSolver/libxcPSolver.a \
	OpenCL/libOpenCL.a \
	bindings/libbindings.a

bigdft_deps =
bigdft_library =

if BUILD_DYNAMIC_LIBS
bigdft_dynamic_library = libbigdft-1.so.@BIGDFT_MINOR_VERSION@.0.@BIGDFT_MICRO_VERSION@
bigdft_soname_library = libbigdft-1.so.@BIGDFT_MINOR_VERSION@
bigdft_deps += @LIB_BIGDFT_DEPS@
bigdft_library += $(bigdft_soname_library) $(bigdft_deps)
AM_LDFLAGS = @FC_RPATH@$(DESTDIR)$(dynlibdir) 
else
bigdft_dynamic_library =
bigdft_soname_library =
bigdft_deps += @LIB_BIGDFT_DEPS@ @LINALG_LIBS@
bigdft_library += libbigdft-1.a $(bigdft_deps)
endif


if BUILD_BINARIES
binaries_sources = bigdft memguess frequencies abscalc utilities
local_binaries = test_forces
scripts_sources = bigdft-tool
sub_dirs_bin = tools tools/bader
EXTRA_PROGRAMS = MDanalysis BigDFT2Wannier WaCo
else
binaries_sources =
local_binaries =
scripts_sources =
sub_dirs_bin =
check_PROGRAMS =
endif

# Additional sources.
# ------------------
if HAVE_LIBCONFIG
sub_dirs_libconfig=convolutions-c
else
sub_dirs_libconfig=
endif

if HAVE_ETSF_IO
etsf_files = wavelib/i-o-etsf.f90
else
etsf_files = wavelib/i-o-etsf_fake.f90
endif


if HAVE_AMBERTOOLS
ambertools_wrap= wrap.c
###for amber supp.:
wrap.o:     wrap.c
	$(NAB) -I$(top_builddir) -c  -g @CC_PIC_FLAG@  $<
else
ambertools_wrap= wrap_fake.f90
endif

##if HAVE_CP2K
##cp2k_interfacemod = modules/module_cp2k.f90
##cp2k_interfacemodCAP = MODULES/MODULE_CP2k.F90
##else
##cp2k_interfacemod = modules/cp2k_fake.f90 modules/module_cp2k.f90
##cp2k_interfacemodCAP = MODULES/CP2K_FAKE.F90 MODULES/MODULE_CP2K_FAKE.F90
##endif

if USE_OCL
OCL_CONV = interface_ocl.f90
else
OCL_CONV =
endif


AM_FCFLAGS = -I. -I$(top_builddir)/includes -I$(srcdir)/wavelib \
	@LIBABINIT_INCLUDE@ @LIB_XC_CFLAGS@ \
	@LIBETSFIO_INCLUDE@ @MPI_INCLUDE@ @LIB_FUTILE_CFLAGS@
AM_CPPFLAGS = @LIB_ARCHIVE_CFLAGS@ -I$(srcdir)/modules -Imodules
#AM_CFLAGS = -Wall -Werror -Wno-strict-aliasing @DGEMMSY_CPPFLAGS@ @CUDA_INCLUDE@ -Wno-deprecated-declarations -Wno-cpp
AM_CFLAGS = @CUDA_INCLUDE@

#objects which depends from module_types
high_level = cluster.f90 \
	init.f90 \
	initialization.f90 \
	input_variables.f90 \
	sumrho.f90 \
	hpsiortho.f90 \
	forces.f90 \
	tail.f90 \
	restart.f90 \
	davidson.f90 \
	ConstrainedDavidson.f90 \
	pdos.f90 \
	tddft.f90 \
	output.f90 \
	external.f90 \
	$(OCL_CONV)

#objects which can be compiled only with module_base
low_level = \
	init/atoms.f90 \
	init/projectors.f90  \
	init/kswfn.f90 \
	init/ionicpot.f90 \
	init/gridmanipulation.f90 \
	init/gauprod.f90 \
	init/denspotd.f90 \
	init/inputguess.f90 \
	init/rcov.f90 \
	init/gautowav.f90 \
	pspconf.f90 \
	init/locreg_init.f90 \
	init/sysprop.f90 \
	init/wavefunctions.f90 \
	init/gaussians_old.f90 \
	init/wfn_init.f90 \
	wfn_opt/precond.f90 \
	wfn_opt/orthogonality.f90 \
	wfn_opt/diis.f90 \
	wfn_opt/applyh.f90 \
	wfn_opt/exctX.f90 \
	wfn_opt/sic.f90 \
	wfn_opt/precond_per_optim.f90 \
	wfn_opt/precond_hyb.f90 \
	wfn_opt/precond_slab.f90 \
	wfn_opt/kernel_per_optim.f90 \
	wfn_opt/kernel_slab_simple.f90 \
	wfn_opt/coupling_matrix.f90 \
	wavelib/transwaves.f90 \
	wavelib/scalprod.f90 \
	wavelib/plotting.f90 \
	wavelib/i-o.f90 \
	wavelib/un-compress.f90 \
	wavelib/gauss_to_daub.f90 \
	$(etsf_files) \
	profiling/memoryestimator.f90 \
	linear/cleanup.f90 \
	linear/filterModule.f90 \
	linear/misc.f90 \
	linear/convolutions_aux.f90 \
	linear/convolutions.f90 \
	linear/locReg_orbitals.f90 \
	linear/linearScaling.f90 \
	linear/initAndUtils.f90 \
	linear/sumrhoLinear.f90 \
	linear/inputguessLinear.f90 \
	linear/orbitalOptimization.f90 \
	linear/mixdiis.f90 \
	linear/copyTypes.f90 \
	linear/derivativeBasis.f90 \
	linear/transferintegrals.f90 \
	linear/cdft.f90 \
	linear/frags.f90\
	$(ambertools_wrap)
##	$(cp2k_interfacemod) 	wfn_opt/exctX_op2p.f90 \

# C objects
c_level = posfiles.c $(DATABASE_SOURCES)

libmain_a_SOURCES = \
	$(c_level) \
	$(low_level) \
	$(high_level)

libabscalc_a_SOURCES = \
	lanczos_base.f90 \
	abscalc_module.f90 \
	lanczos_interface.f90 \
	lanczos.f90 \
	abscalc_init.f90 \
	esatto.f90

pspconf.f90: init/pspconf.in.f90
	$(PYTHON) $(top_srcdir)/config/pspconf.py > psp.inc
	sed '/!!PSP_TABLE!!/r psp.inc' $(srcdir)/init/pspconf.in.f90 > pspconf.f90
	rm -f psp.inc

libbigdft_1_a_SOURCES =

YAML_DATABASE = input_variables_definition.yaml bigdft_bibliography.yaml

#input_variables_description.o: input_variables_definition-inc.h
#build the fortran include file with the stream of input variables into a local string
#the sed command below will transform the entire yaml document into a c string, to be passed to fortran
#input_variables_definition-inc.h: input_variables_definition.yaml
###	$(SED) -e "s/^/\"/;s/$$/\\\n\"/" $< > $@

tmp-libbigdft/extract.stamp: $(subs_static_library) $(out_static_library)
	test -e "tmp-libbigdft" || $(INSTALL) -d -m 755 tmp-libbigdft
	cd tmp-libbigdft ; \
	rm -f *; \
	for lib in $^; do \
	  $(AR) x ../$$lib ; \
	done ; \
	touch extract.stamp

libbigdft-1.a: $(libmain_a_OBJECTS) tmp-libbigdft/extract.stamp
	$(AR) $(ARFLAGS) $@ $(libmain_a_OBJECTS) tmp-libbigdft/*.o
	$(RANLIB) $@

$(bigdft_dynamic_library): $(libmain_a_OBJECTS) tmp-libbigdft/extract.stamp
	$(FC) $(FCFLAGS) -shared $(LDFLAGS) $(AM_LDFLAGS) -Wl,-soname=$(bigdft_soname_library) -o $@ $(libmain_a_OBJECTS) tmp-libbigdft/*.o $(bigdft_deps)
	@chmod a+x $@

libbigdft-1.so.@BIGDFT_MINOR_VERSION@: $(bigdft_dynamic_library)
	ln -fs $^ $@
	ln -fs $@ libbigdft-1.so

install-data-hook:
	if test -n "$(bigdft_dynamic_library)" ; then \
	  cd $(DESTDIR)$(dynlibdir) ; \
	  chmod a+x $(bigdft_dynamic_library) ; \
	  ln -fs $(bigdft_dynamic_library) $(bigdft_soname_library) ; \
	  ln -fs $(bigdft_soname_library) libbigdft-1.so ; \
	fi

#obsolete dependency $(bin_PROGRAMS) $(noinst_PROGRAMS): libbigdft-1.a libdft_common.a

#sources for abscalc exeutable, to be removed from libbigdft.a
#abscalc_additional= 

bigdft_SOURCES = BigDFT.f90
bigdft_LDADD = $(bigdft_library) @FC_EXPORTS@

#rism_SOURCES = rism.f90 rismlowlevel.f90
#rism_LDADD = $(bigdft_library)

#sandbox_SOURCES = sandbox.f90
#sandbox_LDADD = $(bigdft_library)

#oneatom_SOURCES = oneatom.f90
#oneatom_LDADD = $(bigdft_library)

BigDFT2Wannier_SOURCES = BigDFT2Wannier.f90
BigDFT2Wannier_LDADD = $(bigdft_library)

WaCo_SOURCES = WaCo.f90
WaCo_LDADD = $(bigdft_library)

MDanalysis_SOURCES = distances.f90
MDanalysis_LDADD = $(bigdft_library)

abscalc_SOURCES = abscalc.f90 
abscalc_LDADD = libabscalc.a $(bigdft_library)

test_forces_SOURCES = test_forces.f90
test_forces_LDADD = $(bigdft_library)

#to be adjusted: when the compiler has MPI built-in libraries the linker finds multiple definitions
memguess_SOURCES = memguess.f90 MPIfake.f90
memguess_LDADD = $(bigdft_library) @FC_EXPORTS@

utilities_SOURCES = utilities.f90
utilities_LDADD = $(bigdft_library)

frequencies_SOURCES = frequencies.f90
frequencies_LDADD = $(bigdft_library) @FC_EXPORTS@

if BUILD_DYNAMIC_LIBS
bigdft_LDADD += @LIB_FUTILE_LIBS@
memguess_LDADD += @LIB_CHESS_LIBS@ @LIB_PSOLVER_LIBS@ @LIB_FUTILE_LIBS@ @LINALG_LIBS@
frequencies_LDADD += @LIB_FUTILE_LIBS@ @LINALG_LIBS@
abscalc_LDADD += @LIB_PSOLVER_LIBS@ @LINALG_LIBS@
utilities_LDADD += @LIB_CHESS_LIBS@ @LIB_FUTILE_LIBS@
test_forces_LDADD += @LIB_FUTILE_LIBS@
BigDFT2Wannier_LDADD += @LIB_FUTILE_LIBS@ @LINALG_LIBS@
endif

CLEANFILES = *.@MODULE_EXT@ bindings_dbus.* introspect libbigdft-1.*

include $(srcdir)/import_database.mk

clean-local:
	rm -rf tmp-libbigdft

all: all-recursive
	chmod ug+x bigdft-tool
install: install-recursive
	chmod ug+x bigdft-tool

check:
	@if test $(MAKELEVEL) = 0 ; then python $(pythondir)/report.py ; fi

#dependencies
convolutions_aux.o linear/convolutions_aux.o: linear/filterModule.f90 filterModule.o

lanczos_interface.o: abscalc_module.o lanczos_base.o

abscalc_init.o abscalc.o: esatto.o lanczos_interface.o

lanczos.o : lanczos_interface.o

xabsorber.o init/xabsorber.o : init/xabsorber.f90

initAndUtils.o linear/initAndUtils.o: linear/initAndUtils.f90 cleanup.o

linearScaling.o transferintegrals.o: i-o.o

restart.o: i-o.o

# Introspection for Python bindings
# ---------------------------------
if WITH_GOBJECT_INTROSPECTION
BUILT_GIRSOURCES = BigDFT-@BIGDFT_MAJOR_VERSION@.@BIGDFT_MINOR_VERSION@.gir

gir_DATA = $(BUILT_GIRSOURCES)
typelibs_DATA = $(BUILT_GIRSOURCES:.gir=.typelib)

if BUILD_BINDINGS
additonal_sources = $(srcdir)/bindings/bindings.c \
		$(srcdir)/bindings/bindings_regions.c \
		$(srcdir)/bindings/bindings_localfields.c \
		$(srcdir)/bindings/bindings_wavefunctions.c \
		$(srcdir)/bindings/bindings_signals.c \
		$(srcdir)/bindings/bindings_signals_inet.c

endif

$(gir_DATA): $(G_IR_SCANNER) $(bigdft_soname_library)
	echo $$LD_RUN_PATH
	$(G_IR_SCANNER) -v --warn-all \
		--namespace BigDFT \
		--symbol-prefix=bigdft --identifier-prefix=BigDFT_ --identifier-prefix=Bigdft \
		--symbol-prefix=f90 --identifier-prefix=f90_ \
		--nsversion "@BIGDFT_MAJOR_VERSION@.@BIGDFT_MINOR_VERSION@" \
		--add-include-path=$(girdir) \
		--include GObject-2.0 \
		--include GLib-2.0 \
		--include Gio-2.0 \
		--output=$(gir_DATA) \
		--library=bigdft-1 \
		--no-libtool \
		-DGLIB_MAJOR_VERSION=2 -L$(builddir) \
		-I$(srcdir)/bindings -I$(builddir)/bindings \
		$(srcdir)/bindings/bigdft.h \
		$(srcdir)/bindings/bigdft_input_keys.h \
		$(builddir)/bindings/bigdft_cst.h \
		$(srcdir)/bindings/wrappers.c \
		$(srcdir)/bindings/wrapper_inputs.c \
		$(srcdir)/bindings/wrapper_run.c \
		$(srcdir)/bindings/wrapper_atoms.c

#		$(srcdir)/bindings/wrapper_neb.c

$(typelibs_DATA): $(gir_DATA) $(G_IR_COMPILER)
	$(G_IR_COMPILER) \
		--includedir=. \
		--verbose \
		-o $(typelibs_DATA) \
		$(gir_DATA)

CLEANFILES += $(gir_DATA) $(typelibs_DATA)

endif
