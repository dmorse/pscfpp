rpc_tests_solvers_=rpc/tests/solvers/Test.cpp

rpc_tests_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_tests_solvers_:.cpp=.o))

