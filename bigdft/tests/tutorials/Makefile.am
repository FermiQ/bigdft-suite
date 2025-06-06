# Give here the list of existing tests

#SUBDIRS = Xabs/pawpatchgen  Xabs/scfpotgen   Xabs/spectra

EXTRA_TESTDIRS =
SHORT_TESTDIRS =
LONG_TESTDIRS =
if USE_OCL
  EXTRA_TESTDIRS += N2
else
  SHORT_TESTDIRS += N2
endif

SHORT_TESTDIRS += H2O-slab \
	CH4

LONG_TESTDIRS += $(SHORT_TESTDIRS) \
	  CH3- \
	  H2O-bader \
	  scalability \
	  dev

EXTRA_DIST = Figures \
     H1A_first_run.html \
     H1B_grid_cv.html \
     H1C_basic_surfaces.html \
     H1D_wfn_plot.html \
     H2A_MPI_OMP.html \
     H2B_GPU.html \
     H2C-BART.html \
     H2D-splined_saddle.html \
     H2E-minima_hopping.html \
     H2E1-troubleshooting_minima_hopping.html \
     H3B-charge_analysis.html \
     H3C-Wannier.html \
     H3D-Xabs.html \
     H3X-hacking_bigdft.html \
     Mg7-minhop \
     minhop.pdf \
     bart.pdf \
     GPU \
     C2H6 \
     Mg5 \
     Li6 \
     Xabs \
     Wannier


AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
bigdft_library = $(top_builddir)/src/libbigdft-1.so.@BIGDFT_MINOR_VERSION@
else
bigdft_library = $(top_builddir)/src/libbigdft-1.a @LIB_BIGDFT_DEPS@
endif

# Additional rules for the dev tutorial.
dev.check: dev-test/toy_model

dev-test/toy_model: $(abs_top_srcdir)/tests/tutorials/dev/toy_model.f90
	$(FC) -o $@ $(FCFLAGS) $(AM_FCFLAGS) $(abs_top_srcdir)/tests/tutorials/dev/toy_model.f90 \
			$(AM_LDFLAGS) $(LDFLAGS) $(bigdft_library)


toy_model:
	cd ..; $(MAKE) dev/toy_model
dev-N2.out.post-out: dev-N2.out.out toy_model
	rm -f input.perf
	$(run_serial) ./toy_model > log-model.yaml
dev.post-clean:
	rm -f toy_model

# Define the precision for specific directories.
%.report: %.ref
	@case $< in \
          *.memguess.ref | *.out.ref ) mode="--mode=bigdft";; \
          *) mode="";; \
        esac ; \
        case $* in \
          *) prec="1.1e-10" ;; \
        esac ; \
	python @FUTILE_PYTHONDIR@/fldiff.py $$mode --discrepancy=$$prec $*".out" $< | tee $@

%.report.yaml: %.ref.yaml
	@name=`basename $< .out.ref.yaml | sed "s/[^_]*_\?\(.*\)$$/\1/"`  ;\
	tols=`basename $< .out.ref.yaml` ;\
	if test -n "$$name" ; then log=log-$$name.yaml ; else log=log.yaml ; fi ;\
	echo $$log; \
	python @FUTILE_PYTHONDIR@/fldiff_yaml.py \
		-r $< -d $$log -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=$$tols -o $@

include $(srcdir)/../check.mk
