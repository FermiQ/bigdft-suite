AC_PREREQ(2.59)

dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([atlab],[1.0],[luigi.genovese@cea.fr],[atlab])

dnl All temporary building files are put in this directory.
dnl his must be called before AM_INIT_AUTOMAKE
AC_CONFIG_AUX_DIR(config)
dnl Use a config.h file to store macro definitions, required for CUDA.
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_MACRO_DIR([m4])

dnl AM_INIT_AUTOMAKE([-Wno-portability, subdir-objects])
AM_INIT_AUTOMAKE(-Wno-portability)
dnl To keep the object files in the same structure as the source files
dnl AM_INIT_AUTOMAKE([--warnings=no-portability, subdir-objects])

dnl Versioning
ATLAB_MAJOR_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f1)
ATLAB_MINOR_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f2 | cut -d'-' -f1)
ATLAB_MICRO_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f3)
if test -z "$ATLAB_MICRO_VERSION" ; then
  ATLAB_MICRO_VERSION=0
fi
ATLAB_NANO_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f4)
AC_SUBST(ATLAB_MAJOR_VERSION)
AC_SUBST(ATLAB_MINOR_VERSION)
AC_SUBST(ATLAB_MICRO_VERSION)
AC_SUBST(ATLAB_NANO_VERSION)

dnl Check for sed and put it in the SED variable
AC_PROG_SED()

dnl We defined CC for the bindings and CXX for FFTW.
AC_PROG_CC()
AC_PROG_CXX()
AM_PROG_CC_C_O()
AX_CC_RPATH()

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
AX_FC_F2003()
AX_FC_OPENMP()
AX_FC_RPATH()

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

AX_FLIB()

dnl Produce dynamic excecutables.
AX_DYNAMIC_LIBRARIES()

dnl Test Python availability
AX_PYTHON_DEV($ax_build_dynamic)

dnl Test Glib availability
AC_ARG_WITH(gobject, AS_HELP_STRING([--with-gobject], [Compile with GLib (disabled by default).]),
              ax_have_glib=$withval, ax_have_glib="no")
if test x"$ax_have_glib" = x"yes" ; then
  PKG_CHECK_MODULES(GLIB, glib-2.0 gobject-2.0 >= 2.22, [ax_have_glib=yes], [ax_have_glib=no])
fi
if test x"$ax_have_glib" = x"yes" ; then
  AC_DEFINE([HAVE_GLIB], [], [If set, we can call glib.h])
  AC_SUBST([GLIB_PKG], ["glib-2.0"])
  AC_SUBST([GLIB_TRUE], [""])
  AC_SUBST([GLIB_END_TRUE], [""])
  AC_SUBST([GLIB_FALSE], ["/*"])
  AC_SUBST([GLIB_END_FALSE], ["*/"])
else
  AC_SUBST([GLIB_PKG], [""])
  AC_SUBST([GLIB_TRUE], ["/*"])
  AC_SUBST([GLIB_END_TRUE], ["*/"])
  AC_SUBST([GLIB_FALSE], [""])
  AC_SUBST([GLIB_END_FALSE], [""])
fi

AX_ETSF_IO()

AX_OPENBABEL()

dnl Agregate all the optional libs.
LIBS_EXT_DEPS="$LIBETSFIO_LIBS $LIB_FUTILE_LIBS $GLIB_LIBS $PYTHON_LIBS"
AC_SUBST(LIBS_EXT_DEPS, [$LIBS_EXT_DEPS])

dnl if the futile libraries are dynamic add the dependencies
if test $ax_FUTILE_static != "yes" ; then
  LIB_ATLAB_DEPS="$LIB_ATLAB_DEPS $LIBCUDA_LIBS $LINALG_LIBS $MPI_LDFLAGS $MPI_LIBS"
fi
LIB_ATLAB_PKG=""
if test $ax_build_dynamic != "yes" ; then
  LIB_ATLAB_PKG="$LIB_ATLAB_DEPS"
fi
AC_SUBST(LIB_ATLAB_PKG)
AC_SUBST(LIB_ATLAB_DEPS)

LIB_ATLAB_LIBS="-latlab-1"
if test $ax_build_dynamic != "yes" ; then
  LIB_ATLAB_LIBS="$LIB_ATLAB_LIBS $LIBS_EXT_DEPS $LIBS"
fi
AC_SUBST(LIB_ATLAB_LIBS)

AC_SUBST(CONFIG_ARGS, [$ac_configure_args])

dnl Give the name of file.in to transform to file
AC_CONFIG_FILES([
Makefile
atlab.pc
src/Makefile
tests/Makefile
])

AC_OUTPUT

echo "
ATlab configure summary:
========================

Basics:
  Fortran90 compiler:        $FC
  Fortran90 flags:           $FCFLAGS
  Linker flags:              $LDFLAGS
  Linked libraries:          $LIBS_EXT_DEPS $LIBS

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
  With dgemmsy support:      $ax_have_dgemmsy
  With libyaml:              $ax_have_yaml
   | C flags:                $LIB_YAML_CFLAGS
   | linked libraries:       $LIB_YAML_LIBS

Installation paths:
  Source code location:      ${srcdir}
  Prefix:                    $prefix
  Exec prefix:               $exec_prefix
  Binaries:                  $bindir
  Static libraries:          $libdir
  Fortran modules:           $moduledir
  Documentation:             $docdir
"

if test x"$ac_use_mpi" = x"warn" -o x"$withscalapack" = x"warn" ; then
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
