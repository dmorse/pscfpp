rpc_tests_environment_=rpc/tests/environment/Test.cpp

rpc_tests_environment_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_tests_environment_:.cpp=.o))

