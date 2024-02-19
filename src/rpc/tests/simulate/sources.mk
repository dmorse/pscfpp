rpc_tests_simulate_=rpc/tests/simulate/Test.cpp

rpc_tests_simulate_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_tests_simulate_:.cpp=.o))

