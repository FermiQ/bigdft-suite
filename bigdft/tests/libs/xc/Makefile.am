## Process this file with automake to produce Makefile.in

SHORT_CHECK=
LONG_CHECK=xc.out.yaml xc.ref.yaml xc.report.yaml

include $(srcdir)/../check.mk

SUBDIRS =

EXTRA_DIST =

CLEANFILES = xc.out xc.tmp.* xc.report

check_PROGRAMS = test_xc

AM_FCFLAGS = @LIBABINIT_INCLUDE@ -I$(top_builddir)/includes
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src -Wl,-rpath=$(DESTDIR)$(libdir)
bigdft_library = $(top_builddir)/src/libbigdft-1.so.@BIGDFT_MINOR_VERSION@
else
bigdft_library = $(top_builddir)/src/libbigdft-1.a @LIB_BIGDFT_DEPS@
endif

test_xc_SOURCES = test_xc.f90
test_xc_LDADD = $(bigdft_library) @LIB_FUTILE_LIBS@

check: xc.report.yaml

xc.out.yaml: test_xc
	-$(run_parallel) ./test_xc > $@

xc.report.yaml: xc.out.yaml
	python @FUTILE_PYTHONDIR@/fldiff_yaml.py -r $(srcdir)/xc.ref.yaml -d $^ -o $@ \
			-t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=xc
