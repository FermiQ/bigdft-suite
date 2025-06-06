# Give here the list of existing tests
tests_cpu_only = Graphene FeHyb
tests_cpu_only_long = H Mn2 SIC H-Rydberg H-HF
tests_ocl_only = OCL

#initialize
SHORT_TESTDIRS =
LONG_TESTDIRS =
EXTRA_TESTDIRS =

#This will run the automatic tests also for GPU
if USE_CUDA_GPU
  SHORT_TESTDIRS += GPU
else
  EXTRA_TESTDIRS += GPU
endif
if USE_OCL
  SHORT_TESTDIRS += $(tests_ocl_only)
  EXTRA_TESTDIRS += $(tests_cpu_only) $(tests_cpu_only_long)
else
#Add H (box too small) and Graphene (different algorithm) only if we do not use OCL (see README)
#OCL not implemented with spin=4 (Fehyb and Mn2)
  SHORT_TESTDIRS += $(tests_cpu_only)
  LONG_TESTDIRS += $(tests_cpu_only_long)
  EXTRA_TESTDIRS += $(tests_ocl_only)
endif
if HAVE_PSPIO
  SHORT_TESTDIRS += PSPIO
endif

SHORT_TESTDIRS += C \
	taskgroups \
	O2-Spin \
	NLCC \
	PAW-GGA \
        Mg-cubic \
        Na8-multipole \
        C2-chain

LONG_TESTDIRS += $(SHORT_TESTDIRS) \
	  Li+ \
	  H2-spin \
	  SiH4 \
	  SiH4-collocation \
	  NaCl \
	  Ca2 \
	  H2O-CP2K \
	  SiBulk \
	  TiMix \
	  Cavity \
	  ElectField \
	  ElectronGas \
	  CH4-D3 \
	  PAW \
	  counterions

EXTRA_DIST = 

# Give here the pseudo-potentials used by each test.
H.psp: psppar.H
Li+.psp: psppar.Li
H2-spin.psp: HGH/psppar.H
SiH4.psp: HGH-K/psppar.H HGH-K/psppar.Si
TiMix.psp: HGH-K/psppar.Ti HGH/psppar.Na
TiSmear.psp: HGH-K/psppar.Ti HGH/psppar.Na
NaCl.psp: HGH/psppar.Na HGH/psppar.Cl
Ca2.psp: psppar.Ca
Mn2.psp: psppar.Mn
H2O-CP2K.psp: HGH/psppar.H HGH/psppar.O
SiBulk.psp: HGH-K/psppar.Si
AgBulk.psp: HGH/psppar.Ag
Graphene.psp: psppar.C
FeHyb.psp: HGH/psppar.Fe
GPU.psp: HGH/psppar.Si
OCL.psp: HGH/psppar.Si
O2-Spin.psp: HGH-K/psppar.O
ElectField.psp: HGH-K/psppar.H HGH-K/psppar.O
PAW-GGA.psp: C.GGA_PBE-JTH.xml
PSPIO.psp: Li.UPF Si.pbe-hgh.UPF

SIC-test.post-in:
	if [ ! -d data-PZ ] ; then mkdir data-PZ; fi
	mv wavefunctions.gau data-PZ/

# Define the precision for specific directories.
%.report: %.ref
	@case $< in \
          *.memguess.ref | *.out.ref) mode="--mode=bigdft";; \
          *) mode="";; \
        esac ; \
        case $* in \
          *)                    prec="1.1e-10" ;; \
        esac ; \
    echo "python @FUTILE_PYTHONDIR@/fldiff.py $$mode --discrepancy=$$prec $*.out $<" ; \
	python @FUTILE_PYTHONDIR@/fldiff.py $$mode --discrepancy=$$prec $*".out" $< | tee $@

%.report.yaml: %.ref.yaml
	@name=`basename $< .out.ref.yaml | sed "s/[^_]*_\?\(.*\)$$/\1/"`; \
	tols=`basename $< .out.ref.yaml` ;\
	if test -n "$$name" ; then log=log-$$name.yaml ; else log=log.yaml ; fi ;\
	echo "python @FUTILE_PYTHONDIR@/fldiff_yaml.py -r $< -d $$log -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=$$tols -o $@" ; \
	python @FUTILE_PYTHONDIR@/fldiff_yaml.py \
		-r $< -d $$log -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=$$tols -o $@ 

H-Rydberg-test.post-in:
	@if test -f "H.perf" ; then \
		echo >> H.yaml ; \
		echo "perf:" >> H.yaml ; \
		awk '{ printf("  %s: %s\n",$$1,$$2); }' H.perf >> H.yaml ; \
		echo "  multipole_preserving: Yes" >> H.yaml ; \
		rm H.perf; \
	fi

include $(srcdir)/../../check.mk


