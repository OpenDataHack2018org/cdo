# Process this file with autoconf to produce a configure script.

AC_INIT(cdo, 1.0.4, Uwe.Schulzweida@zmaw.de)

CONFIG_ABORT=yes

AC_ARG_VAR(AR,     [AR archive command]) 

AC_CONFIG_AUX_DIR(config)
AC_CANONICAL_TARGET

AM_CONFIG_HEADER(src/config.h)
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CHECK_PROGS(AR, ar, [ar])

# Check large file support on 32 bit systems
AC_SYS_LARGEFILE

# Checks for library functions.
AC_FUNC_MMAP

# Checks for library malloc.
AC_CHECK_LIB(malloc, mallinfo)

# Checks for the availability of ANSI-C99 functions
AC_CHECK_FUNCS(isnan)

# Checks for compiler
COMPILER="$CC $CFLAGS"
AC_DEFINE_UNQUOTED(COMPILER, ["$COMPILER"], [Compiler])

# Check compiler version
case "$CC" in
  gcc)    COMP_VERSION=`$CC --version | head -n 1`;;
  g++)    COMP_VERSION=`$CC --version | head -n 1`;;
  sxc++)  COMP_VERSION=`$CC -V 2>&1   | tail -n 1`;;
  xlc*)   COMP_VERSION=`$CC    2>&1   | head -n 1`;;
  *)      COMP_VERSION=`$CC -V 2>&1   | head -n 1`;;
esac

if test -z "$COMP_VERSION" ; then COMP_VERSION="unknown"; fi;
AC_DEFINE_UNQUOTED(COMP_VERSION, ["$COMP_VERSION"], [Compiler version])

# Checks for username and hostname
USERNAME=$LOGNAME
if test -z "$USERNAME" ; then USERNAME=$USER; fi;
if test -z "$USERNAME" ; then USERNAME="unknown"; fi;

AC_DEFINE_UNQUOTED(USER_NAME, ["$USERNAME"], [User name])

if test -z "$HOST" ; then HOST="unknown"; fi;

AC_DEFINE_UNQUOTED(HOST_NAME, ["$HOST"], [Host name])

#  ----------------------------------------------------------------------
#  Checks pthread library
#
AC_ARG_WITH([threads],
            [AC_HELP_STRING([--with-threads],
                            [Use the pthread library [default=yes]])],,)

if test "$with_threads" = "no" ; then
  AC_MSG_CHECKING([for pthread library])
  AC_MSG_RESULT([suppressed])
else
  AC_CHECK_HEADERS(pthread.h)
  AC_CHECK_LIB(pthread, pthread_create)
fi

#  ----------------------------------------------------------------------
#  Compile interface with internal GRIB library
#
AC_ARG_WITH([grib],
            [AC_HELP_STRING([--with-grib],
                            [Use the internal GRIB library [default=yes]])],,
            AC_DEFINE(HAVE_LIBGRIB, [1], [Define to 1 for GRIB interface]))

if test -z "$with_grib" ; then with_grib=yes ; fi

AC_MSG_CHECKING([for GRIB support])
AC_MSG_RESULT($with_grib)

#  ----------------------------------------------------------------------
#  Compile interface with internal SERVICE library
#
AC_ARG_WITH([service],
            [AC_HELP_STRING([--with-service],
                            [Use the internal SERVICE library [default=yes]])],,
            AC_DEFINE(HAVE_LIBSERVICE, [1], [Define to 1 for SERVICE interface]))

if test -z "$with_service" ; then with_service=yes ; fi

AC_MSG_CHECKING([for SERVICE support])
AC_MSG_RESULT($with_service)

#  ----------------------------------------------------------------------
#  Compile interface with internal EXTRA library
#
AC_ARG_WITH([extra],
            [AC_HELP_STRING([--with-extra],
                            [Use the internal EXTRA library [default=yes]])],,
            AC_DEFINE(HAVE_LIBEXTRA, [1], [Define to 1 for EXTRA interface]))

if test -z "$with_extra" ; then with_extra=yes ; fi

AC_MSG_CHECKING([for EXTRA support])
AC_MSG_RESULT($with_extra)

#  ----------------------------------------------------------------------
#  Compile interface with internal IEG library
#
AC_ARG_WITH([ieg],
            [AC_HELP_STRING([--with-ieg],
                            [Use the internal IEG library [default=yes]])],,
            AC_DEFINE(HAVE_LIBIEG, [1], [Define to 1 for IEG interface]))

if test -z "$with_ieg" ; then with_ieg=yes ; fi

AC_MSG_CHECKING([for IEG support])
AC_MSG_RESULT($with_ieg)

#  --with-xxx from HDF configure
#  ----------------------------------------------------------------------
#  Fake --with-xxx option to allow us to create a help message for the
#  following --with-xxx options which can take either a =DIR or =INC,LIB
#  specifier.
# 
AC_ARG_WITH([fnord],
  [
 For the following --with-xxx options, you can specify where the header
 files and libraries are in two different ways:

    --with-xxx=INC,LIB - Specify individually the include directory and
                         library directory separated by a comma
    --with-xxx=DIR     - Specify only the directory which contains the
                         include/ and lib/ subdirectories
  ])

#  ----------------------------------------------------------------------
#  Is the hdf5 present? It has a header file `hdf5.h' and a library
#  `-lhdf5' and their locations might be specified with the `--with-hdf5'
#  command-line switch. The value is an include path and/or a library path.
#  If the library path is specified then it must be preceded by a comma.
#
AC_ARG_WITH([hdf5],
            [AC_HELP_STRING([--with-hdf5=DIR],
                            [Use the HDF5 library [default=no]])],,withval=no)

case $withval in
  yes)
    HAVE_HDF5="yes"
    AC_CHECK_HEADERS(hdf5.h)
    AC_CHECK_LIB(hdf5, H5Fopen,, unset HAVE_HDF5)

    if test -z "$HAVE_HDF5" -a -n "$CONFIG_ABORT"; then
      AC_MSG_ERROR([couldn't find hdf5 library])
    fi
    ;;
  no)
    HAVE_HDF5="no"
    AC_MSG_CHECKING([for hdf5 library])
    AC_MSG_RESULT([suppressed])
    ;;
  *)
    HAVE_HDF5="yes"
    case "$withval" in
      *,*)
        hdf5_inc="`echo $withval |cut -f1 -d,`"
        hdf5_lib="`echo $withval |cut -f2 -d, -s`"
        ;;
      *)
        if test -n "$withval"; then
          hdf5_inc="$withval/include"
          hdf5_lib="$withval/lib"
        fi
        ;;
    esac

#   Trying to include -I/usr/include and -L/usr/lib is redundant and
#   can mess some compilers up.
    if test "X$hdf5_inc" = "X/usr/include"; then
      hdf5_inc=""
    fi
    if test "X$hdf5_lib" = "X/usr/lib"; then
      hdf5_lib=""
    fi

#   Trying to include -I/usr/local/include and -L/usr/local/lib is redundant and
#   can mess some compilers up.
    if test "X$hdf5_inc" = "X/usr/local/include"; then
      hdf5_inc=""
    fi
    if test "X$hdf5_lib" = "X/usr/local/lib"; then
      hdf5_lib=""
    fi

    saved_CPPFLAGS="$CPPFLAGS"
    saved_LDFLAGS="$LDFLAGS"

    if test -n "$hdf5_inc"; then
      CPPFLAGS="$CPPFLAGS -I$hdf5_inc"
    fi

    AC_CHECK_HEADERS(hdf5.h,, CPPFLAGS="$saved_CPPFLAGS")

    if test -n "$hdf5_lib"; then
      LDFLAGS="$LDFLAGS -L$hdf5_lib"
    fi

    AC_CHECK_LIB(hdf5, H5Fopen,, LDFLAGS="$saved_LDFLAGS"; unset HAVE_HDF5)

    if test -z "$HAVE_HDF5" -a -n "$CONFIG_ABORT"; then
      AC_MSG_ERROR(couldn't find hdf5 library)
    fi
    ;;
esac

#  ----------------------------------------------------------------------
# Begin DAP (aka DODS) support
# Store pre-DAP LIBS value to restore from later
LIBS_no_DAP="${LIBS}"

AC_ARG_ENABLE(dap,AS_HELP_STRING([--enable-dap],[Build OPeNDAP/DODS-enabled clients if possible [[default=no]]]),enable_dap=${enableval},enable_dap=no)

if test "${enable_dap}" = "yes"; then
   if which ncdap-config > /dev/null; then
      CPPFLAGS="${CPPFLAGS} `ncdap-config --cflags`"
      nc_dap_ver=`ncdap-config --version`
      LIBS="`ncdap-config --libs` ${LIBS}"	
   else
      AC_MSG_ERROR(couldn't find ncdap-config)     
   fi
fi

#  ----------------------------------------------------------------------
#  Is the netcdf present? It has a header file `netcdf.h' and a library
#  `-lnetcdf' and their locations might be specified with the `--with-netcdf'
#  command-line switch. The value is an include path and/or a library path.
#  If the library path is specified then it must be preceded by a comma.
#
AC_ARG_WITH([netcdf],
            [AC_HELP_STRING([--with-netcdf=DIR],
                            [Use the netCDF library [default=no]])],,withval=no)

case $withval in
  yes)
    HAVE_NETCDF="yes"
    AC_CHECK_HEADERS(netcdf.h)
    AC_CHECK_LIB(netcdf, nc_open,, unset HAVE_NETCDF)

    if test -z "$HAVE_NETCDF" -a -n "$CONFIG_ABORT"; then
      AC_MSG_ERROR([couldn't find netcdf library])
    fi
    ;;
  no)
    HAVE_NETCDF="no"
    AC_MSG_CHECKING([for netcdf library])
    AC_MSG_RESULT([suppressed])
    ;;
  *)
    HAVE_NETCDF="yes"
    case "$withval" in
      *,*)
        netcdf_inc="`echo $withval |cut -f1 -d,`"
        netcdf_lib="`echo $withval |cut -f2 -d, -s`"
        ;;
      *)
        if test -n "$withval"; then
          netcdf_inc="$withval/include"
          netcdf_lib="$withval/lib"
        fi
        ;;
    esac

#   Trying to include -I/usr/include and -L/usr/lib is redundant and
#   can mess some compilers up.
    if test "X$netcdf_inc" = "X/usr/include"; then
      netcdf_inc=""
    fi
    if test "X$netcdf_lib" = "X/usr/lib"; then
      netcdf_lib=""
    fi

#   Trying to include -I/usr/local/include and -L/usr/local/lib is redundant and
#   can mess some compilers up.
    if test "X$netcdf_inc" = "X/usr/local/include"; then
      netcdf_inc=""
    fi
    if test "X$netcdf_lib" = "X/usr/local/lib"; then
      netcdf_lib=""
    fi

    saved_CPPFLAGS="$CPPFLAGS"
    saved_LDFLAGS="$LDFLAGS"

    if test -n "$netcdf_inc"; then
      CPPFLAGS="$CPPFLAGS -I$netcdf_inc"
    fi

    AC_CHECK_HEADERS(netcdf.h,, CPPFLAGS="$saved_CPPFLAGS")

    if test -n "$netcdf_lib"; then
      LDFLAGS="$LDFLAGS -L$netcdf_lib"
    fi

    if test "${enable_dap}" = "yes"; then
      AC_CHECK_LIB(nc-dap, nc_open,, LDFLAGS="$saved_LDFLAGS"; unset HAVE_NETCDF)
if test $ac_cv_lib_nc_dap_nc_open = yes; then
cat >>confdefs.h <<_ACEOF
#define HAVE_LIBNETCDF 1
_ACEOF
fi
    else
      AC_CHECK_LIB(netcdf, nc_open,, LDFLAGS="$saved_LDFLAGS"; unset HAVE_NETCDF)
    fi

    if test -z "$HAVE_NETCDF" -a -n "$CONFIG_ABORT"; then
      AC_MSG_ERROR(couldn't find netcdf library)
    fi
    ;;
esac

AC_SUBST([CPPFLAGS])
AC_SUBST([FC],$F77)
AC_SUBST([AR])

AC_OUTPUT(Makefile src/Makefile)
