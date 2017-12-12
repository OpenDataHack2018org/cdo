
noinst_PROGRAMS = cdo_process_creation_test operator_1_to_1 operator_1_to_2 chained_operators
# bin_PROGRAMS = cdo_process_creation_test operator_1_to_1 operator_1_to_2


tests_CPPFLAGS = -I$(top_srcdir)/test/bandit_tests -I$(top_srcdir)/test/bandit_tests/bandit
tests_LDADD = $(top_builddir)/src/libcdo.la $(top_builddir)/libcdi/src/libcdi.la

cdo_process_creation_test_SOURCES = process_creation_tests.cc
cdo_process_creation_test_CPPFLAGS = $(tests_CPPFLAGS)
cdo_process_creation_test_LDADD = $(tests_LDADD)

operator_1_to_1_SOURCES = operator_1_to_1.bandit.cc
operator_1_to_1_CPPFLAGS = $(tests_CPPFLAGS)
operator_1_to_1_LDADD = $(tests_LDADD)

operator_1_to_2_SOURCES = operator_1_to_2.bandit.cc
operator_1_to_2_CPPFLAGS = $(tests_CPPFLAGS)
operator_1_to_2_LDADD = $(tests_LDADD)

chained_operators_SOURCES = chained_operators.cc
chained_operators_CPPFLAGS = $(tests_CPPFLAGS)
chained_operators_LDADD = $(tests_LDADD)
LD_LIBRARY_PATH = $(tests_LDADD)
TESTS = ./cdo_process_creation_test ./operator_1_to_1 ./operator_1_to_2 ./chained_operators

