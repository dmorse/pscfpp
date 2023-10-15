pscf_tests_=pscf/tests/Test.cc

pscf_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_:.cc=.o))

