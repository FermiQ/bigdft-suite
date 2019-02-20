SUBDIRS = xcPSolver \
	OP2P \
	psp \
	IO-wf \
	wrappers \
	overlapgeneral \
	coordinatetransform \
	xc

if USE_BLACS
SUBDIRS += dsygv
endif 

if USE_CUDA_GPU
SUBDIRS += CUDA
endif

if USE_OCL
SUBDIRS += OpenCL
endif

if BUILD_BINDINGS
SUBDIRS += ftypes
endif

check_PROGRAMS = PS_Check

AM_FCFLAGS = -I$(top_builddir)/includes -I. @LIBABINIT_INCLUDE@ @LIB_XC_CFLAGS@  @MPI_INCLUDE@
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src -Wl,-rpath=$(DESTDIR)$(libdir)
bigdft_library = $(top_builddir)/src/libbigdft-1.so.@BIGDFT_MINOR_VERSION@
else
bigdft_library = $(top_builddir)/src/libbigdft-1.a @LIB_BIGDFT_DEPS@
endif

PS_Check_SOURCES = PS_Check.f90
PS_Check_LDADD = $(bigdft_library) @LIB_PSOLVER_LIBS@ @LIB_FUTILE_LIBS@ @LINALG_LIBS@

CLEANFILES = *.out.yaml *.report.yaml report_remarks *dat fort* log.yaml time.yaml

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report: 
	python $(pythondir)/report.py

regtest: PS_Check
	python $(pythondir)/f_regtest.py --fldiff=$(pythondir)/fldiff_yaml.py \
	       --tols=$(abs_srcdir)/tols.yaml --srcdir=$(srcdir) \
	       --exclusive=${CHECK_ONLY}

define F_REGTEST_INSTRUCTIONS
- PS_Check:
   runs:
    - $(run_parallel) ./PS_Check --ndim=[57,48,63] --ixc=0  --geocode=F  > PS_Check.out.yaml
    - $(run_parallel) ./PS_Check --ndim=[57,48,63] --ixc=1  --geocode=F >> PS_Check.out.yaml
    - $(run_parallel) ./PS_Check --ndim=[57,48,63] --ixc=11 --geocode=F >> PS_Check.out.yaml
    - $(run_parallel) ./PS_Check --ndim=[57,48,63] --ixc=13 --geocode=F >> PS_Check.out.yaml
    - $(run_parallel) ./PS_Check --ndim=[64,64,64] --ixc=0  --geocode=P >> PS_Check.out.yaml
    - $(run_parallel) ./PS_Check --ndim=[64,64,64] --ixc=1  --geocode=P >> PS_Check.out.yaml
    - $(run_parallel) ./PS_Check --ndim=[64,64,64] --ixc=11 --geocode=P >> PS_Check.out.yaml
    - $(run_parallel) ./PS_Check --ndim=[64,64,64] --ixc=13 --geocode=P >> PS_Check.out.yaml
    - $(run_parallel) ./PS_Check --ndim=[32,64,48] --ixc=0  --geocode=S >> PS_Check.out.yaml
    - $(run_parallel) ./PS_Check --ndim=[32,64,48] --ixc=1  --geocode=S >> PS_Check.out.yaml
    - $(run_parallel) ./PS_Check --ndim=[32,64,48] --ixc=11 --geocode=S >> PS_Check.out.yaml
    - $(run_parallel) ./PS_Check --ndim=[32,64,48] --ixc=13 --geocode=S >> PS_Check.out.yaml
    - cp -f PS_Check.out.yaml tmp
    - cat tmp | grep -v 'Unable to read mpd.hosts' > PS_Check.out.yaml
    - rm -f tmp

endef
export F_REGTEST_INSTRUCTIONS
