## Process this file with automake to produce Makefile.in

SUBDIRS = \
	  cubic \
	  postSCF 

if USE_OCL
else
  SUBDIRS += linear
endif

EXTRA_DIST = \
	README \
	inguess.dat

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report: 
	python @FUTILE_PYTHONDIR@/report.py
