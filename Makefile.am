## Process this file with automake to produce Makefile.in
SUBDIRS =
if USE_INTERNAL_CDI
SUBDIRS += libcdi
endif
SUBDIRS += src contrib
#
EXTRA_DIST=config/default OPERATORS doc/cdo.pdf doc/cdo_refcard.pdf cdo.spec
#
ACLOCAL_AMFLAGS = -I m4
#
CLEANFILES  = `ls *~`
