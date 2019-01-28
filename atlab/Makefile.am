
EXTRA_DIST= m4

ACLOCAL_AMFLAGS = -I m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = atlab.pc

SUBDIRS = src tests

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report:
	python $(top_builddir)/tests/report.py

clean-local:
	rm -rf includes
