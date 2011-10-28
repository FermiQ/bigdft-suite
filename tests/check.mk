# Generic part, of the testing Makefiles.
# Possible calls are:
#  make in: generate all input dirs.
#  make failed-check: run check again on all directories with missing report
#                     or failed report.
#  make X.in: generate input dir for directory X.
#  make X.check: generate a report for directory X (if not already existing).
#  make X.recheck: force the creation of the report in directory X.
#  make X.clean: clean the given directroy X.

if USE_MPI
  mpirun_message=mpirun
else
  mpirun_message=
endif
if USE_OCL
oclrun_message = oclrun
else
oclrun_message =
endif

AM_FCFLAGS = -I$(top_builddir)/src -I$(top_builddir)/src/PSolver -I$(top_builddir)/src/modules @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@

PSPS = psppar.H \
       psppar.C \
       psppar.Li \
       psppar.Ca \
       psppar.Mn \
       psppar.N \
       psppar.Si \
       HGH/psppar.H \
       HGH/psppar.Na \
       HGH/psppar.Cl \
       HGH/psppar.O \
       HGH/psppar.Si \
       HGH/psppar.Fe \
       HGH/psppar.Mg \
       HGH/psppar.Ag \
       HGH/psppar.N \
       HGH-K/psppar.H \
       HGH-K/psppar.Si \
       HGH-K/psppar.N \
       HGH-K/psppar.O \
       HGH-K/psppar.Ti \
       Xabs/psppar.Fe

INS = $(TESTDIRS:=.in)
RUNS = $(TESTDIRS:=.run)
CHECKS = $(TESTDIRS:=.check)
FAILEDCHECKS = $(TESTDIRS:=.recheck)
CLEANS = $(TESTDIRS:=.clean)

in: $(INS)

check: $(CHECKS) report

clean: $(CLEANS)

distclean : $(CLEANS)
	rm -rf Makefile

failed-check: $(FAILEDCHECKS) report

report:
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

%.memguess.out: $(abs_top_builddir)/src/memguess
	name=`basename $@ .out.out | sed "s/[^_]*_\?\(.*\)$$/\1/"` ; \
	$(abs_top_builddir)/src/memguess 1 $$name > $@
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.out.out: $(abs_top_builddir)/src/bigdft
	name=`basename $@ .out.out | sed "s/[^_]*_\?\(.*\)$$/\1/"` ; \
	$(run_parallel) $(abs_top_builddir)/src/bigdft $$name > $@
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.geopt.mon.out:
	$(MAKE) -f ../Makefile $*.out.out && mv geopt.mon $@
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.dipole.dat.out: %.out.out
	$(abs_top_builddir)/src/tools/bader/bader data/electronic_density.cube > bader.out && mv dipole.dat $@
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.freq.out: $(abs_top_builddir)/src/frequencies
	$(run_parallel) $(abs_top_builddir)/src/frequencies > $@
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.NEB.out: $(abs_top_builddir)/src/NEB NEB_include.sh NEB_driver.sh
	rm -f triH.NEB.it*
	$(abs_top_builddir)/src/NEB < input | tee $@
	rm -rf triH.NEB.0*
	rm -f gen_output_file velocities_file
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.splsad.out: $(abs_top_builddir)/src/splsad
	$(run_parallel) $(abs_top_builddir)/src/splsad > $@
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.minhop.out: $(abs_top_builddir)/src/global
	$(run_parallel) $(abs_top_builddir)/src/global > $@
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.xabs.out: $(abs_top_builddir)/src/abscalc
	name=`basename $@ .xabs.out | sed "s/[^_]*_\?\(.*\)$$/\1/"` ; \
	$(abs_top_builddir)/src/abscalc $$name > $@
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"

$(PSPS):
	ln -fs $(abs_top_srcdir)/utils/PSPfiles/$@ 

%.clean:
	@dir=`basename $@ .clean` ; \
	rm -f $$dir.* ; \
        if test x"$(srcdir)" = x"." ; then \
	  cd $$dir ; \
	  for i in psppar.* ; do \
	    if test -L $i ; then \
	      rm -f $i ; \
	    fi ; \
	  done ; \
          rm -f *.out *.mon *.report default* *.prc; \
	  rm -fr data data-*; \
	  rm -f anchorpoints* fort.* nogt.* path*.xyz vogt.* ; \
	  rm -f latest.pos.force.*.dat fort.* CPUlimit test ; \
	  rm -f cheb_spectra_* alphabeta* b2B_xanes.* local_potentialb2B* ; \
	  $(MAKE) -f ../Makefile $$dir".post-clean"; \
        else \
          rm -rf $$dir ; \
        fi ; \
        echo "Clean in "$$dir

%.post-in: ;
%.psp: ;
%.post-clean: ;
%.post-out: ;

NEB_driver.sh:
	ln -fs $(abs_top_srcdir)/src/$@

in_message:
	@if test -n "$(run_ocl)" ; then \
	  echo "==============================================" ; \
	  echo "Will generate a 'input.perf' file to force OCL" ; \
	  echo "==============================================" ; \
	fi

$(INS): in_message
	@dir=`basename $@ .in` ; \
        if ! test x"$(srcdir)" = x"." ; then \
          if [ ! -d $$dir ] ; then mkdir $$dir ; fi ; \
          for i in $(srcdir)/$$dir/* ; do cp -f $$i $$dir; done ; \
        fi ; \
	if ! test -f $(srcdir)/$$dir/input.perf ; then \
	  rm -f $$dir/input.perf ; \
	  if test -n "$(run_ocl)" ; then \
	    echo "ACCEL OCLGPU" > $$dir/input.perf ; \
	  fi ; \
	fi ; \
        cd $$dir && $(MAKE) -f ../Makefile $$dir".psp"; \
        $(MAKE) -f ../Makefile $$dir".post-in"; \
        echo "Input prepared in "$$dir" dir. make $$dir.run available"
	touch $@

run_message:
	@if test -n "$(run_parallel)" ; then \
	  echo "==============================================" ; \
	  echo "Will run tests in parallel with '$$run_parallel'" ; \
	  echo "==============================================" ; \
	fi

%.run: %.in run_message
	@dir=`basename $@ .run` ; \
        runs="$(srcdir)/$$dir/*.ref" ; \
	tgts=`for r in $$runs ; do echo $$(basename $$r .ref)".out"; done` ; \
        cd $$dir && $(MAKE) -f ../Makefile $$tgts ; \
        echo "Tests have run in "$$dir" dir. make $$dir.check available"
	touch $@

%.check: %.run
	@dir=`basename $@ .check` ; \
        chks="$(srcdir)/$$dir/*.ref" ; \
	tgts=`for c in $$chks ; do echo $$(basename $$c .ref)".report"; done` ; \
        cd $$dir && $(MAKE) -f ../Makefile $$tgts
	touch $@

%.recheck: %.in
	@dir=`basename $@ .recheck` ; \
        refs="$$dir/*.ref" ; \
	for r in $$refs ; do \
	  rep=`basename $$r .ref`".report" ; \
	  if ! grep -qs "succeeded\|passed" $$dir/$$rep ; then \
	    target=` basename $$r .ref` ; \
	    rm -f $$dir/$$target".out" $$dir/$$target".report" ; \
	    cd $$dir && $(MAKE) -f ../Makefile $$target".out" $$target".report" && cd - ; \
	  fi \
	done
	touch $*".check"

# Avoid copying in dist the builddir files.
distdir: $(DISTFILES)
	@srcdirstrip=`echo "$(srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	topsrcdirstrip=`echo "$(top_srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	list='$(DISTFILES)'; \
	  dist_files=`for file in $$list; do echo $$file; done | \
	  sed -e "s|^$$srcdirstrip/||;t" \
	      -e "s|^$$topsrcdirstrip/|$(top_builddir)/|;t"`; \
	case $$dist_files in \
	  */*) $(MKDIR_P) `echo "$$dist_files" | \
			   sed '/\//!d;s|^|$(distdir)/|;s,/[^/]*$$,,' | \
			   sort -u` ;; \
	esac; \
	for file in $$dist_files; do \
	  d=$(srcdir); \
	  if test -d $$d/$$file; then \
	    dir=`echo "/$$file" | sed -e 's,/[^/]*$$,,'`; \
	    if test -d "$(distdir)/$$file"; then \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    if test -d $(srcdir)/$$file && test $$d != $(srcdir); then \
	      cp -fpR $(srcdir)/$$file "$(distdir)$$dir" || exit 1; \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    cp -fpR $$d/$$file "$(distdir)$$dir" || exit 1; \
	  else \
	    test -f "$(distdir)/$$file" \
	    || cp -p $$d/$$file "$(distdir)/$$file" \
	    || exit 1; \
	  fi; \
	done

# Doc messages.
all:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) foot_message ; fi

head_message:
	@echo "========================================================="
	@echo " This is a directory for tests. Beside the 'make check'"
	@echo " one can use the following commands:"
	@echo "  make in:           generate all input dirs."
	@echo "  make failed-check: run check again on all directories"
	@echo "                     with missing report or failed report."
	@echo "  make X.in:         generate input dir for directory X."
	@echo "  make X.check:      generate a report for directory X"
	@echo "                     (if not already existing)."
	@echo "  make X.recheck:    force the creation of the report in"
	@echo "                     directory X."
	@echo "  make X.clean:      clean the given directroy X."

mpirun: head_message
	@echo ""
	@echo " Use the environment variable run_parallel"
	@echo "     ex: export run_parallel='mpirun -np 2'  "

oclrun: head_message $(mpirun_message)
	@echo ""
	@echo " Use the environment variable run_ocl"
	@echo "     ex: export run_ocl='on' to use OpenCL acceleration"

foot_message: $(mpirun_message) $(oclrun_message) head_message
	@echo "========================================================="

