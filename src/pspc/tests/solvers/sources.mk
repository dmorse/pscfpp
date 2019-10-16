pspc_tests_solvers_=pspc/tests/solvers/Test.cc

pspc_tests_solvers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspc_tests_solvers_))
pspc_tests_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_tests_solvers_:.cc=.o))

