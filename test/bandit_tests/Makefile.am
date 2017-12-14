
noinst_PROGRAMS = chained_operators variable_inputs_operators
# bin_PROGRAMS = cdo_process_creation_test operator_1_to_1 operator_1_to_2


tests_CPPFLAGS = -I$(top_srcdir)/test/bandit_tests -I$(top_srcdir)/test/bandit_tests/bandit
tests_LDADD = $(top_builddir)/src/libcdo.la $(top_builddir)/libcdi/src/libcdi.la

chained_operators_SOURCES = chained_operators.cc
chained_operators_CPPFLAGS = $(tests_CPPFLAGS)
chained_operators_LDADD = $(tests_LDADD)

variable_inputs_operators_SOURCES = variable_inputs_operators.cc
variable_inputs_operators_CPPFLAGS = $(tests_CPPFLAGS)
variable_inputs_operators_LDADD = $(tests_LDADD)
LD_LIBRARY_PATH = $(tests_LDADD)
TESTS = ./chained_operators ./variable_inputs_operators

