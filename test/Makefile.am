CLEANFILES =

export

# COLOR ALWAYS also if no color capable terminal detected
#     this is needed i.e. to run the tests via a oneline ssh cmd:
#     ssh host '... make check'
AM_COLOR_TESTS=always

TEST_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(top_srcdir)/config/tap-driver.sh


# tests which should pass
TESTS = threads.test tsformat.test wildcard.test File.test Comp.test Compc.test \
        Read_grib.test Read_netcdf.test Copy_netcdf.test Cat.test Gridarea.test  \
        MapReduce.test Ninfo.test Gradsdes.test  \
        Detrend.test Genweights.test Remap.test Remap2.test Remap3.test Remapeta.test \
        Select.test Spectral.test Ensstat.test Enspctl.test Gridboxstat.test \
        Vertstat.test Fldstat.test Fldpctl.test Merstat.test Zonstat.test \
        Timstat.test Timstat2.test Timselstat.test Seasstat.test Timpctl.test \
        Runstat.test Runpctl.test Multiyearstat.test Ydrunstat.test \
        Mergetime.test EOF.test \
        Vertint.test Afterburner.test Arithc.test Arith.test Ymonarith.test Expr.test \
        Collgrid.test Setmiss.test Smooth.test \
        Filter.test

# tests which should fail
XFAIL_TESTS = 

#        $(top_srcdir)/test/test_Remap.sh \
#        $(top_srcdir)/test/test_info.py
#        $(top_srcdir)/test/test_diff.py
#        $(top_srcdir)/test/test_Arith.py
#        $(top_srcdir)/test/test_Arithc.py
#        $(top_srcdir)/test/test_Selvar.py
#        $(top_srcdir)/test/test_intgridbil.py

#EXTRA_DIST = $(TESTS) $(top_srcdir)/test/testStreams.py

CDO        = $(top_builddir)/src/cdo

DATAPATH   = $(top_srcdir)/test/data

PYTHONPATH = $(top_srcdir)/contrib/python:$(top_srcdir)/test

#.PHONY: $(TESTS)

CLEANFILES += `ls *.pyc *.debug cerr* cout* ifile* ofile*`

#AUTOMAKE_OPTIONS = color-tests
