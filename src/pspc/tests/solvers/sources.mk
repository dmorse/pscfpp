pspc_tests_solvers_=pspc/tests/solvers/Test.cpp

pspc_tests_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_tests_solvers_:.cpp=.o))

