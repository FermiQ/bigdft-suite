AC_PREREQ(2.59)

dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([Chebyshev Sparse Solvers],[0.2.4],[stephan.mohr@bsc.es],[CheSS])

AC_CONFIG_MACRO_DIR([m4])

dnl All temporary building files are put in this directory.
dnl his must be called before AM_INIT_AUTOMAKE
AC_CONFIG_AUX_DIR(config)
dnl Use a config.h file to store macro definitions, required for CUDA.
AC_CONFIG_HEADERS(config.h)

dnl Init something, don't know exactly what...
AM_INIT_AUTOMAKE([-Wno-portability])

dnl Versioning
CHESS_MAJOR_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f1)
CHESS_MINOR_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f2 | cut -d'-' -f1)
CHESS_MICRO_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f3)
if test -z "$CHESS_MICRO_VERSION" ; then
  CHESS_MICRO_VERSION=0
fi
CHESS_NANO_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f4)
AC_SUBST(CHESS_MAJOR_VERSION)
AC_SUBST(CHESS_MINOR_VERSION)
AC_SUBST(CHESS_MICRO_VERSION)
AC_SUBST(CHESS_NANO_VERSION)

dnl Check for sed and put it in the SED variable
AC_PROG_SED()

dnl We defined CC for the bindings.
AC_PROG_CC()
AM_PROG_CC_C_O()
dnl If FC environment variable is not already set, it
dnl looks for a modern Fortran compiler, prefering newest
dnl fortran available (e. g. Fortran95). Then FC variable
dnl is set.
AC_PROG_FC([], [Fortran 90])

dnl Defining Fortran language for tests.
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)

AC_FC_WRAPPERS()
AX_FC_DEFS()
AX_FC_OPENMP()

dnl Set if the module name are capitalized or not with the FC compiler.
AX_FC_MOD()
AM_CONDITIONAL(CAPITALIZE, test "$ax_fc_mod_capitalize" = "yes")
AC_SUBST(MODULE_EXT, $ax_fc_mod_ext)
dnl default installation directories
AX_FC_MODULEDIR()

dnl Set the FCFLAGS and FFLAGS variable
if test -z "$FFLAGS" ; then
  AC_SUBST(FFLAGS, "-O2")
fi
if test -n "$FCFLAGS_f90" ; then
  FCFLAGS="$FCFLAGS $FCFLAGS_f90"
fi
if test -z "$FCFLAGS" ; then
  AC_SUBST(FCFLAGS, "-O2")
fi

dnl Look for ranlib and canonicalize it with the $host variable
AC_CHECK_TOOL(RANLIB, ranlib, :)
if test "$ranlib" = ":" ; then
  AC_MSG_ERROR(["No 'ranlib' program found."])
fi
dnl Look for ar and canonicalize it with the $host variable
AC_CHECK_TOOL(AR, ar, :)
if test "$ar" = ":" ; then
  AC_MSG_ERROR(["No 'ar' program found."])
fi

dnl Enable or not MPI.
AX_MPI()

dnl Python YAML support.
AX_PYYAML()

dnl The utilities executable
AC_ARG_ENABLE(binaries, AS_HELP_STRING([--disable-binaries], [build the binaries (default = yes).]), ac_build_bin_utilities=$enableval, ac_build_bin_utilities="yes")
AM_CONDITIONAL(BUILD_BINARIES, test "$ac_build_bin_utilities" = "yes")
AC_MSG_CHECKING([for CheSS binaries build])
AC_MSG_RESULT([$ac_build_bin_utilities])


AX_ATLAB()

dnl Produce dynamic executables.
AX_DYNAMIC_LIBRARIES([], [f_lib_initialize], [$LIB_FUTILE_LIBS])

dnl Test accelerators
AX_ACC_CUDA()

dnl PEXSI
AC_ARG_WITH(pexsi, AS_HELP_STRING([--with-pexsi],
            [Give the link-line for PEXSI.]),
            ac_pexsi_linkline=$withval, ac_pexsi_linkline=)
ax_have_pexsi="no"
if test -n "$ac_pexsi_linkline" ; then
  PEXSI_LINKLINE="$ac_pexsi_linkline"
  ax_have_pexsi="yes"
fi
AM_CONDITIONAL(HAVE_PEXSI, test $ax_have_pexsi = "yes")

dnl MINPACK
AC_ARG_WITH(minpack, AS_HELP_STRING([--with-minpack],
            [Give the link-line for MINPACK.]),
            ac_minpack_linkline=$withval, ac_minpack_linkline=)
ax_have_minpack="no"
if test -n "$ac_minpack_linkline" ; then
  MINPACK_LINKLINE="$ac_minpack_linkline"
  ax_have_minpack="yes"
fi
AM_CONDITIONAL(HAVE_MINPACK, test $ax_have_minpack = "yes")

AX_LINALG()

dnl Agregate all the optional libs.
LIB_CHESS_DEPS="$LIB_ATLAB_LIBS $PEXSI_LINKLINE"
dnl if the futile libraries are dynamic add the dependencies
if test $ax_FUTILE_static != "yes" ; then
  LIB_CHESS_DEPS="$LIB_CHESS_DEPS $LIBCUDA_LIBS $LINALG_LIBS $MPI_LDFLAGS $MPI_LIBS"
fi
LIB_CHESS_PKG=""
if test $ax_build_dynamic != "yes" ; then
  LIB_CHESS_PKG="$LIB_CHESS_DEPS"
fi
AC_SUBST(LIB_CHESS_PKG)
AC_SUBST(LIB_CHESS_DEPS)

dnl Test for the developers' documentation.
AX_DOC()

AC_SUBST(CONFIG_ARGS, [$ac_configure_args])

dnl Give the name of file.in to transform to file
AC_CONFIG_FILES([
Makefile
chess.pc
chess_environment.sh
chessvars.sh
doc/Makefile
doc/Doxyfile
sparsematrix/Makefile
src/configure.inc
src/Makefile
src/python/Makefile
src/python/CheSS/Makefile
chesspexsi/Makefile
tests/Makefile
])

AC_OUTPUT

echo "
Chebyshev Sparse Solvers configure summary:
================================

Basics:
  Fortran90 compiler:        $FC
  Fortran90 flags:           $FCFLAGS
  Fortran77 flags:           $FFLAGS
  Linker flags:              $LDFLAGS
  Linked libraries:          $LIB_CHESS_DEPS $LIBS

Build:
  Dynamic library:           $ax_build_dynamic

Options:
  With OpenMP:               $ax_fc_openmp_msg
  With MPI:                  $ac_use_mpi
   | include dir.:           $MPI_INCLUDE
   | linker flags:           $MPI_LDFLAGS
   | linked libraries:       $MPI_LIBS
   | MPI_INIT_THREAD avail:  $ac_with_mpinithread
   | MPI2 support:           $ac_with_mpi2
  With PEXSI:                $ax_have_pexsi
  With MINPACK:              $ax_have_minpack

Installation paths:
  Source code location:      ${srcdir}
  Prefix:                    $prefix
  Exec prefix:               $exec_prefix
  Binaries:                  $bindir
  Static libraries:          $libdir
  Fortran modules:           $moduledir
  Documentation:             $docdir
"

if test x"$ac_use_mpi" = x"warn" -o x"$ax_have_etsf_io" = x"warn" -o x"$withscalapack" = x"warn" ; then
  echo "Warnings:"
  if test x"$ac_use_mpi" = x"warn" ; then
    echo "  MPI:                       disabled by tests"
  fi
  if test x"$withscalapack" = x"warn" ; then
    echo "  Scalapack binding:         disabled by tests"
  fi

fi

dnl Do not touch this line. On the Cray machine in Basel, the preprocessor
dnl in Fortran sources does not allow C comments (/*  */), so we remove them.
$SED -n '1h;1!H;${;g;s%/\*[[^*]]*\*/%%g;p;}' config.h > config.inc
