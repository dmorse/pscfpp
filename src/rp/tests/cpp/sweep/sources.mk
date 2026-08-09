rpc_tests_sweep_=rpc/tests/sweep/Test.cpp

rpc_tests_sweep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_tests_sweep_:.cpp=.o))

