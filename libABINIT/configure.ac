

AC_PREREQ(2.59)

dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([ABINIT for BigDFT - a subpart of ABINIT for BigDFT],[6.8+12],[Damien.Caliste@cea.fr],[libABINIT])

dnl Init something, don't know exactly what...
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE(-Wno-portability)

dnl Check for sed and put it in the SED variable
AC_PROG_SED()
dnl We defined CC for the ABINIT part.
AC_PROG_CC()

dnl If FC environment variable is not already set, it
dnl looks for a modern Fortran compiler, prefering newest
dnl fortran available (e. g. Fortran95). Then FC variable
dnl is set.
AC_PROG_FC([], [Fortran 90])
AC_FC_SRCEXT(f90)

AX_FC_DEFS()

dnl Set if the module name are capitalized or not with the FC compiler.
AX_FC_MOD()
AM_CONDITIONAL(CAPITALIZE, test "$ax_fc_mod_capitalize" = "yes")
AC_SUBST(MODULE_EXT, $ax_fc_mod_ext)
dnl default installation directories
AX_FC_MODULEDIR()

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

AX_MPI()
AX_DYNAMIC_LIBRARIES([], [f_lib_initialize], [$LIB_FUTILE_LIBS])
dnl Test the libPAW library.
ac_internal_libpaw="yes"
AC_ARG_WITH(internal-libpaw, AS_HELP_STRING([--with-internal-libpaw],
                                            [compile also libPAW (default = yes).]),
            ac_internal_libpaw=$withval, ac_internal_libpaw="yes")
AM_CONDITIONAL(BUILD_LIBPAW, test x"$ac_internal_libpaw" = x"yes")

AC_ARG_WITH(libpaw-host, AS_HELP_STRING([--with-libpaw-host],
                                        [choice of low-level for LibPAW (default = BigDFT).]),
            ac_libpaw_host=$withval, ac_libpaw_host="BigDFT")

LIBPAW_HOST_DEFS=""
if test x"$ac_internal_libpaw" = x"yes" ; then
  dnl FLib detection for libPAW.
  AX_FLIB()
  dnl libXC detection for libPAW.
  AX_LIBXC()
  dnl Choose the libPAW host mode.
  AC_MSG_CHECKING([for libPAW host macro])
  case "$ac_libpaw_host" in
    bigdft | BigDFT) LIBPAW_HOST_DEFS="-DHAVE_LIBPAW_BIGDFT";;
    abinit | ABINIT) LIBPAW_HOST_DEFS="-DHAVE_LIBPAW_ABINIT";;
    *) LIBPAW_HOST_DEFS="";;
  esac
  AC_MSG_RESULT([$LIBPAW_HOST_DEFS])
fi
AC_SUBST(LIBPAW_HOST_DEFS)

moduledir=${includedir}
AC_SUBST(moduledir)

dnl Give the name of file.in to transform to file
AC_CONFIG_FILES([
Makefile
src/Makefile
])

AC_OUTPUT
