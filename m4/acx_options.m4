AC_DEFUN([ACX_CDO_OPTIONS],
[
#  ----------------------------------------------------------------------
#  Checks for multithreaded compiling + linking
ENABLE_THREADS=no
AC_ARG_WITH([threads],
            [AC_HELP_STRING([--with-threads=<yes/no/directory>],
                            [Compile + link for multithreading [default=yes]])],
            [],
            [with_threads=yes])
THREADS_INCLUDE=''
THREADS_LIBS=''
AS_CASE([$with_threads],
        [no],[AC_MSG_CHECKING([multithreading])
              AC_MSG_RESULT([suppressed])],
        [yes],[AX_PTHREAD([AC_DEFINE([HAVE_LIBPTHREAD],[1],[Define 1 for multithread support])
                           ENABLE_THREADS=yes],[AC_MSG_ERROR([multithreaded settings NOT found])])
               LIBS="$PTHREAD_LIBS $LIBS"
               CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
               CXXFLAGS="$CXXFLAGS $PTHREAD_CFLAGS"
               CC="$PTHREAD_CC"
               AS_ECHO(["CC:$CC CFLAGS:$CFLAGS LIBS:$LIBS"])],
        [*],[THREADS_ROOT=$with_threads
             LDFLAGS="-L$THREADS_ROOT/lib $LDFLAGS"
             CPPFLAGS="-I$THREADS_ROOT/include $CPPFLAGS "
             AC_CHECK_HEADERS(pthread.h)
             AC_CHECK_LIB([pthread],[pthread_create])
             ENABLE_THREADS=yes
             THREADS_LIBS=" -L$THREADS_ROOT/lib -lpthread"
             THREADS_INCLUDE=" -I$THREADS_ROOT/include"])
AC_SUBST([ENABLE_THREADS])
AC_SUBST([THREADS_INCLUDE])
AC_SUBST([THREADS_LIBS])
#  ----------------------------------------------------------------------
#  Link application with HDF5 library, required for netcdf4
HDF5_ROOT=''
HDF5_INCLUDE=''
HDF5_LIBS=''
AC_ARG_WITH([hdf5],
            [AS_HELP_STRING([--with-hdf5=<yes|no|directory> (default=no)],[location of HDF5 library])],
            [AS_CASE(["$with_hdf5"],
                     [no],[AC_MSG_CHECKING([for hdf5 library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([hdf5.h])
                            AC_SEARCH_LIBS([H5Fopen],
                                           [hdf5],
                                           [AC_DEFINE([HAVE_LIBHDF5],[1],[Define to 1 for HDF5 support])],
                                           [AC_MSG_ERROR([Cannot link to hdf5 library!])])
                            AC_SEARCH_LIBS([H5DSis_scale],
                                           [hdf5_hl],
                                           [have_hdf5_hl=yes],
                                           [AC_MSG_NOTICE([Cannot find hdf5 high level interface!])
                                            have_hdf5_hl=no])
                            AS_IF([test "x$have_libhdf5_hl" = xyes],
                                  [HDF5_LIBS=" -lhdf5_hl -lhdf5"],
                                  [HDF5_LIBS=" -lhdf5"])
                            ],
                     [*],[AS_IF([test -d "$with_hdf5"],
                                [HDF5_ROOT="$with_hdf5"
                                 LDFLAGS="-L$HDF5_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$HDF5_ROOT/include $CPPFLAGS"
                                 AC_CHECK_HEADERS([hdf5.h])
                                 AC_SEARCH_LIBS([H5Fopen],
                                                [hdf5],
                                                [AC_DEFINE([HAVE_LIBHDF5],[1],[Define to 1 for HDF5 support])],
                                                [AC_MSG_ERROR([Cannot link to hdf5!])])
                                 AC_SEARCH_LIBS([H5DSis_scale],
                                                [hdf5_hl],
                                                [have_hdf5_hl=yes],
                                                [AC_MSG_NOTICE([Cannot link to hdf5 high level interface!])
                                                have_hdf5_hl=no])
                                 AS_IF([test "x$have_libhdf5_hl" = 'xyes'],
                                       [HDF5_LIBS=" -L$HDF5_ROOT/lib -lhdf5_hl -lhdf5"],
                                       [HDF5_LIBS=" -L$HDF5_ROOT/lib -lhdf5"])
                                 HDF5_INCLUDE=" -I$HDF5_ROOT/include"],
                                [AC_MSG_NOTICE([$with_hdf5 is not a directory! HDF5 suppressed])])])],
            [AC_MSG_CHECKING([for hdf5 library])
             AC_MSG_RESULT([suppressed])])
AC_SUBST([HDF5_ROOT])
AC_SUBST([HDF5_INCLUDE])
AC_SUBST([HDF5_LIBS])

#  ----------------------------------------------------------------------
#  Link application with UDUNITS2 library
AC_ARG_WITH([udunits2],
            [AS_HELP_STRING([--with-udunits2=<directory>],
                            [Specify location of UDUNITS2 library.])],
            [AS_CASE(["$with_udunits2"],
                     [no],[AC_MSG_CHECKING([for udunits2 library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([udunits2.h])
                            AC_SEARCH_LIBS([ut_parse],
                                           [udunits2],
                                           [AC_DEFINE([HAVE_LIBUDUNITS2],[1],[Define to 1 for UDUNITS2 support])],
                                           [AC_MSG_ERROR([Could not link to udunits2 library!])])
                            AC_SUBST([UDUNITS_LDFLAGS],[" -ludunits2"])
                            AC_SUBST([UDUNITS_INCLUDE],[""])],
                     [*],[UDUNITS_ROOT=$with_udunits2
                          AS_IF([test -d "$UDUNITS_ROOT"],
                                [LDFLAGS="$LDFLAGS -L$UDUNITS_ROOT/lib"
                                 CPPFLAGS="$CPPFLAGS -I$UDUNITS_ROOT/include -I$UDUNITS_ROOT/include/udunits2"
                                 AC_CHECK_HEADERS([udunits2.h])
                                 AC_CHECK_HEADERS([udunits2/udunits2.h])
                                 AC_SEARCH_LIBS([ut_parse],
                                                [udunits2],
                                                [AC_DEFINE([HAVE_LIBUDUNITS2],[1],[Define to 1 for UDUNITS2 support])],
                                                [AC_MSG_ERROR([Could not link to udunits2 library!])])
                                 AC_SUBST([UDUNITS_LDFLAGS],[" -L$UDUNITS_ROOT/lib -ludunits2"])
                                 AC_SUBST([UDUNITS_INCLUDE],[" -I$UDUNITS_ROOT/include"])],
                                [AC_MSG_ERROR([$UDUNITS_ROOT is not a directory! UDUNITS2 suppressed])])])],
            [AC_MSG_CHECKING([for the UDUNITS2 library])
             AC_MSG_RESULT([suppressed])])

# -----------------------------------------------------------------------
# test for UUID libraries needed by CMOR
# (util-linux libuuid in 3.1.1 and before, OSSP UUID in 3.1.2 and after)
ACX_UUID
#  ----------------------------------------------------------------------
#  Link application with CMOR library
CMOR_LIBS=''
AC_ARG_WITH([cmor],
            [AS_HELP_STRING([--with-cmor=<directory>],[Specify location of CMOR library.])],
            [AS_CASE(["$with_cmor"],
                     [no],[AC_MSG_CHECKING([for cmor library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([cmor.h])
                            LIBS="${LIBS+${LIBS} }$UUID_C_LIB"
                            AC_SEARCH_LIBS([cmor_load_table],[cmor],[AC_DEFINE([HAVE_LIBCMOR],[1],[Define to 1 for CMOR support])],
                                           [AC_MSG_ERROR([Could not link to cmor library!])])
                            CMOR_LIBS="-lcmor $UUID_C_LIB"],
                     [*],[CMOR_ROOT=$with_cmor
                          AS_IF([test -d "$CMOR_ROOT"],
                                [LDFLAGS="$LDFLAGS -L$CMOR_ROOT/lib"
                                 LIBS="${LIBS+${LIBS} }$UUID_C_LIB"
                                 CPPFLAGS="$CPPFLAGS -I$CMOR_ROOT/include -I$CMOR_ROOT/include/cdTime -I$CMOR_ROOT/include/json-c"
                                 AC_SEARCH_LIBS([cmor_load_table],
                                                [cmor],
                                                [AC_DEFINE([HAVE_LIBCMOR],[1],[Define to 1 for CMOR support])],
                                                [AC_MSG_ERROR([Could not link to cmor library!])])
                                 CMOR_LIBS="-L$CMOR_ROOT/lib -lcmor $UUID_C_LIB"],
                                [AC_MSG_ERROR([$CMOR_ROOT is not a directory! CMOR suppressed])])])],
            [AC_MSG_CHECKING([for the CMOR library])
             AC_MSG_RESULT([suppressed])])
AC_SUBST([CMOR_LIBS])
#  ----------------------------------------------------------------------
#  Compile with fftw support
AC_MSG_CHECKING([for FFTW3 support])
AC_ARG_WITH([fftw3],
    [AS_HELP_STRING([--with-fftw3],
      [enable support for fftw3])],
    [],
    [with_fftw3=no])

  AS_IF([test "x$with_fftw3" != xno],
      [AC_CHECK_HEADERS([fftw3.h])
    AC_SEARCH_LIBS([fftw_cleanup],[fftw3],
      [AC_DEFINE([HAVE_LIBFFTW3],[1],[FFTW3 library is present if defined to 1])],
      [AC_MSG_RESULT([Could not link to fftw3 library])])])

#  ----------------------------------------------------------------------
#  Checks for PROJ.4 library
AC_ARG_WITH([proj],
            [AS_HELP_STRING([--with-proj=<directory>],
                            [Specify location of PROJ library for cartographic projections.])],
            [AS_CASE(["$with_proj"],
                     [no],[AC_MSG_CHECKING([for proj library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([proj_api.h])
                            AC_SEARCH_LIBS([pj_init],[proj],[AC_DEFINE([HAVE_LIBPROJ],[1],[Define to 1 for PROJ support])],
                                           [AC_MSG_ERROR([Could not link to PROJ library!])])
                            AC_SUBST([PROJ_LDFLAGS],[" -lproj"])
                            AC_SUBST([PROJ_INCLUDE],[""])],
                     [*],[PROJ_ROOT=$with_proj
                          AS_IF([test -d "$PROJ_ROOT"],
                                [LDFLAGS="-L$PROJ_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$PROJ_ROOT/include $CPPFLAGS"
                                 AC_CHECK_HEADERS([proj_api.h])
                                 AC_SEARCH_LIBS([pj_init],
                                                [proj],
                                                [AC_DEFINE([HAVE_LIBPROJ],[1],[Define to 1 for PROJ support])],
                                                [AC_MSG_ERROR([Could not link to PROJ library!])])
                                 AC_SUBST([PROJ_LDFLAGS],[" -L$PROJ_ROOT/lib -lproj"])
                                 AC_SUBST([PROJ_INCLUDE],[" -I$PROJ_ROOT/include"])],
                                [AC_MSG_ERROR([$PROJ_ROOT is not a directory! PROJ suppressed])])])],
            [AC_MSG_CHECKING([for the PROJ library])
             AC_MSG_RESULT([suppressed])])
#  ----------------------------------------------------------------------
#  Checks for CURL library
AC_ARG_WITH([curl],
            [AS_HELP_STRING([--with-curl=<directory>],
                            [Specify location of CURL library.])],
            [AS_CASE(["$with_curl"],
                     [no],[AC_MSG_CHECKING([for curl library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([curl/curl.h])
                            AC_SEARCH_LIBS([curl_global_init],[curl],[AC_DEFINE([HAVE_LIBCURL],[1],[Define to 1 for CURL support])],
                                           [AC_MSG_ERROR([Could not link to CURL library!])])
                            AC_SUBST([CURL_LDFLAGS],[" -lcurl"])
                            AC_SUBST([CURL_INCLUDE],[""])],
                     [*],[CURL_ROOT=$with_curl
                          AS_IF([test -d "$CURL_ROOT"],
                                [LDFLAGS="-L$CURL_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$CURL_ROOT/include $CPPFLAGS"
                                 AC_CHECK_HEADERS([curl/curl.h])
                                 AC_SEARCH_LIBS([curl_global_init],
                                                [curl],
                                                [AC_DEFINE([HAVE_LIBCURL],[1],[Define to 1 for CURL support])],
                                                [AC_MSG_ERROR([Could not link to CURL library!])])
                                 AC_SUBST([CURL_LDFLAGS],[" -L$CURL_ROOT/lib -lcurl"])
                                 AC_SUBST([CURL_INCLUDE],[" -I$CURL_ROOT/include"])],
                                [AC_MSG_ERROR([$CURL_ROOT is not a directory! CURL suppressed])])])],
            [AC_MSG_CHECKING([for the CURL library])
             AC_MSG_RESULT([suppressed])])
#  ----------------------------------------------------------------------
#  Compile application with MAGICS (xml required)
MAGICS_ROOT=''
MAGICS_INCLUDE=''
MAGICS_LIBS=''
AC_ARG_WITH([magics],
            [AS_HELP_STRING([--with-magics=<yes|no|directory>],[location of magics library (lib and include subdirs)])],
            [AS_CASE(["$with_magics"],
                     [no],[AC_MSG_CHECKING([for magics library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([magics_api.h])
                            AC_SEARCH_LIBS([mag_open],
                                           [MagPlus],
                                           [AC_DEFINE([HAVE_LIBMAGICS],[1],[Define to 1 for MAGICS support])],
                                           [AC_MSG_ERROR([Could not link to magics library])])
                            AC_SUBST([MAGICS_LIBS],[" -lMagPlus"])],
                     [*],[AS_IF([test -d "$with_magics"],
                                [MAGICS_ROOT=$with_magics
                                 LDFLAGS="-L$MAGICS_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$MAGICS_ROOT/include/magics $CPPFLAGS"
                                 AC_CHECK_HEADERS([magics_api.h])
                                 AC_SEARCH_LIBS([mag_open],
                                                [MagPlus],
                                                [AC_DEFINE([HAVE_LIBMAGICS],[1],[Define to 1 for MAGICS support])],
                                                [AC_MSG_ERROR([Could not link to magics library])])
                                 MAGICS_LIBS=" -L$MAGICS_ROOT/lib -lMagPlus"
                                 MAGICS_INCLUDE=" -I$MAGICS_ROOT/include/magics"],
                                [AC_MSG_NOTICE([$with_magics is not a directory! MAGICS suppressed])])])],
            [AC_MSG_CHECKING([for MAGICS library])
             AC_MSG_RESULT([suppressed])])
AC_SUBST([MAGICS_ROOT])
AC_SUBST([MAGICS_INCLUDE])
AC_SUBST([MAGICS_LIBS])
AC_ARG_WITH([libxml2],
            [AS_HELP_STRING([--with-libxml2=<yes|no|directory>],[location of libxml2 library (lib and include subdirs)])],
            [AS_CASE(["$with_libxml2"],
                     [no],[AC_MSG_CHECKING([for libxml2 library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[CPPFLAGS="$CPPFLAGS -I/usr/include/libxml2"
                            AC_CHECK_HEADERS([libxml/parser.h])
                            AC_CHECK_HEADERS([libxml/tree.h])
                            AC_SEARCH_LIBS([xmlInitMemory],
                                           [xml2],
                                           [AC_DEFINE([HAVE_LIBXML2],[1],[Define to 1 for XML2 support])],
                                           [AC_MSG_ERROR([Could not link to libxml2 library])])
                            AC_SUBST([XML2_LIBS],[" -lxml2"])],
                     [*],[AS_IF([test -d "$with_libxml2"],
                                [XML2_ROOT=$with_libxml2
                                 LDFLAGS="-L$XML2_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$XML2_ROOT/include/libxml2 $CPPFLAGS"
                                 AC_CHECK_HEADERS([libxml/parser.h])
                                 AC_CHECK_HEADERS([libxml/tree.h])
                                 AC_SEARCH_LIBS([xmlInitMemory],
                                                [xml2],
                                                [AC_DEFINE([HAVE_LIBXML2],[1],[Define to 1 for XML2 support])],
                                                [AC_MSG_ERROR([Could not link to libxml2 library])])
                                 XML2_LIBS=" -L$XML2_ROOT/lib -lxml2"
                                 XML2_INCLUDE=" -I$XML2_ROOT/include/libxml2"],
                                [AC_MSG_NOTICE([$with_libxml2 is not a directory! XML2 suppressed])])])],
            [AC_MSG_CHECKING([for XML2 library])
             AC_MSG_RESULT([suppressed])])
AM_CONDITIONAL([ENABLE_MAGICS],[test ! x$with_magics = 'x' -a ! x$with_libxml2 = 'x'])
#  ----------------------------------------------------------------------
#  How to build CDI into CDO? 
INTERNAL_CDI_DIR=libcdi
# At the moment, there are two possible CDI bindings
# (default)          linking directly to CDI's object files, i.e. a libtool
#                    convenience library
# (--enable-cdi-lib) build and link to a shared CDI library
AC_MSG_CHECKING([for build a separate CDI library and link CDO to it])
AC_ARG_ENABLE([cdi-lib],
              [AS_HELP_STRING([--enable-cdi-lib],[build, install and link to a CDI library [default=no]])],
              [AS_IF([test "x$enable_cdi_lib" != "xno"],
                     [enable_cdi_lib=yes],
                     [enable_cdi_lib=no
                      CDO_DISABLE_CDILIB=1
                      export CDO_DISABLE_CDILIB])],
              [enable_cdi_lib=no
               CDO_DISABLE_CDILIB=1
               export CDO_DISABLE_CDILIB])
AC_MSG_RESULT([$enable_cdi_lib])
# save CDI binding mode for later automake use
AM_CONDITIONAL([ENABLE_CDI_LIB],[test x$enable_cdi_lib = 'xyes'])
# create shell variables for the representation of configure results
AS_IF([test x$enable_cdi_lib = 'xno'],[AC_SUBST([ENABLE_CDI_LIB],[false])],[AC_SUBST([ENABLE_CDI_LIB],[true])])
# scan for CDI as a subproject
AC_CONFIG_SUBDIRS([libcdi])
#  ----------------------------------------------------------------------
#  Build a static CDO
AC_MSG_CHECKING([for building an additional static CDO binary])
AC_ARG_ENABLE([all-static],
              [AS_HELP_STRING([--enable-all-static],[build a completely statically linked CDO binary [default=no]])],
              [AS_IF([test "x$enable_all_static" != "xno"],
                     [enable_all_static=yes],
                     [enable_all_static=no])],
              [enable_all_static=no])
AC_MSG_RESULT([$enable_all_static])
AM_CONDITIONAL([ENABLE_ALL_STATIC],[test x$enable_all_static = 'xyes'])
#  ----------------------------------------------------------------------
#  Build CDO with HIRLAM extensions
AC_MSG_CHECKING([for HIRLAM extensions])
AC_ARG_ENABLE([hirlam-extensions],
              [AS_HELP_STRING([--enable-hirlam-extensions],[HIRLAM extensions [default=no]])],
              [AS_IF([test "x$enable_hirlam_extensions" != "xno"],
                    [AC_DEFINE(HIRLAM_EXTENSIONS,[1],[Define to 1 for HIRLAM extensions])
                     enable_hirlam_extensions=yes],
                    [enable_hirlam_extensions=no])],
              [enable_hirlam_extensions=no])
AC_MSG_RESULT([$enable_hirlam_extensions])
AM_CONDITIONAL([ENABLE_HIRLAM_EXTENSIONS],[test x$enable_hirlam_extensions = 'xyes'])
#
])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
