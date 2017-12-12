# Process this file with automake to produce Makefile.in
SUBDIRS = libcdi src/lib/ncl src contrib test/data test
if MAINTAINER_MODE
SUBDIRS += test/bandit_tests
endif
#
EXTRA_DIST=config/default OPERATORS doc/cdo.pdf doc/cdo_cmor.pdf doc/cdo_eca.pdf doc/cdo_magics.pdf doc/cdo_refcard.pdf cdo.spec README
#
ACLOCAL_AMFLAGS = -I m4
#
CLEANFILES  = `ls *~ 2> /dev/null`
