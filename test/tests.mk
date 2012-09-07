TESTFILES := $(wildcard *.py)
TESTS     := $(patsubst %.py,%,$(TESTFILES))

define testTarget_template
$(1): $(top_srcdir)/test/$(1).py
	echo CDO="$(top_builddir)/src/cdo" PYTHONPATH="$(top_srcdir)/contrib/python" python $(top_srcdir)/test/$(1).py > $(top_builddir)/test/$(1); chmod u+x $(top_builddir)/test/$(1); ./$(top_builddir)/test/$(1)
endef
$(foreach test,$(TESTS),$(eval $(call testTarget_template,$(test))))


listTests:
	@echo $(TESTS)

test: $(TESTS)

.PHONY: $(TESTS)
CLEANFILES += $(TESTS)
