
EXTRA_DIST= m4

ACLOCAL_AMFLAGS = -I m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = atlab.pc

SUBDIRS = src tests

DISTCLEANFILES = config.inc

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report:
	python @FUTILE_PYTHONDIR@/report.py

clean-local:
	rm -rf includes
