# Process this file with automake to produce Makefile.in
SUBDIRS = libcdi src/lib/ncl src/lib/yac src contrib test/data test
if MAINTAINER_MODE
SUBDIRS += test/bandit_tests
endif
#
ACLOCAL_AMFLAGS = -I m4 -I libcdi/m4
#
EXTRA_DIST=config/default OPERATORS doc/cdo.pdf doc/cdo_cmor.pdf doc/cdo_eca.pdf doc/cdo_magics.pdf doc/cdo_refcard.pdf cdo.spec README
#
#
CLEANFILES  = `ls *~ 2> /dev/null`
