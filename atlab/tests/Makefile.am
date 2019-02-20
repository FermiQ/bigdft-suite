EXTRA_DIST = numerics_check.ref.yaml babel.ref.yaml domain.yaml TiO2.cif

AM_LDFLAGS = -L$(top_builddir)/src

if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += @FC_RPATH@$(abs_top_builddir)/src @FC_RPATH@$(DESTDIR)$(libdir)
endif

if HAVE_OPENBABEL
babel_prg = babel
babel_report = babel.out.report.yaml
if BUILD_DYNAMIC_LIBS
else
AM_LDFLAGS += @FC_RPATH@$(abs_top_builddir)/src
endif
endif

check_PROGRAMS = numerics_check $(babel_prg)

check: regtest report

report:
	@if test $(MAKELEVEL) = 0 ; then python $(pythondir)/report.py ; fi

CLEANFILES = *.out.yaml *.report.yaml report_remarks report outfile.xyz OpenBabel.err

numerics_check_SOURCES = numerics_check.f90
numerics_check_LDADD = @LIB_ATLAB_LIBS@ @LIB_FUTILE_LIBS@

babel_SOURCES = babel.f90
babel_LDADD = @LIB_ATLAB_LIBS@ @LIB_FUTILE_LIBS@

AM_FCFLAGS = -I. @LIB_FUTILE_CFLAGS@ @LIBETSFIO_INCLUDE@

numerics_check.o: $(top_builddir)/src/libatlab-1.a

regtest: numerics_check $(babel_prg)
	python $(pythondir)/f_regtest.py --fldiff=$(pythondir)/fldiff_yaml.py --tols=$(srcdir)/tols.yaml --srcdir=$(srcdir)

define F_REGTEST_INSTRUCTIONS
- numerics_check:
   runs:
    - cp -f $(srcdir)/domain.yaml $(builddir)
    - ./numerics_check
if HAVE_OPENBABEL
- babel:
   runs:
    - ./babel -i $(srcdir)/TiO2.cif
endif
endef
export F_REGTEST_INSTRUCTIONS
