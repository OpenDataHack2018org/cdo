CLEANFILES =

export

TESTS = $(top_srcdir)/test/test_info.py \
	$(top_srcdir)/test/test_diff.py \
	$(top_srcdir)/test/test_intgridbil.py

EXTRA_DIST = $(TESTS)

CDO   = $(top_builddir)/src/cdo

PYTHONPATH = $(top_srcdir)/contrib/python

listTests:
	@echo $(TESTS)

test: $(TESTS)

.PHONY: $(TESTS)

CLEANFILES += `ls *.pyc`
