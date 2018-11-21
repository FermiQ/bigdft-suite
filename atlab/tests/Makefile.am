EXTRA_DIST = \
	numerics_check.ref.yaml

check_PROGRAMS = numerics_check

check: regtest report

report:
  @if test $(MAKELEVEL) = 0 ; then python $(top_builddir)/tests/report.py ; fi

CLEANFILES = *.out.yaml *.report.yaml *.@MODULE_EXT@ report_remarks

numerics_check_SOURCES = numerics_check.f90
numerics_check_LDADD = @LIB_FUTILE_LIBS@

numerics_check.o: $(top_builddir)/src/libatlab-1.a

regtest: numerics_check
	python $(pythondir)/f_regtest.py --fldiff=$(pythondir)/fldiff_yaml.py --tols=$(srcdir)/tols.yaml --srcdir=$(srcdir)

define F_REGTEST_INSTRUCTIONS
- numerics_check:
  runs:
  - ./numerics_check
endef
export F_REGTEST_INSTRUCTIONS
