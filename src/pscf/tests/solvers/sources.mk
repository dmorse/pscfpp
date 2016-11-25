pscf_tests_solvers_=pscf/tests/solvers/Test.cc

pscf_tests_solvers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pscf_tests_solvers_))
pscf_tests_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_solvers_:.cc=.o))

