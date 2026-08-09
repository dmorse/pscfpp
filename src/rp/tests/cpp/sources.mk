rpc_tests_=rpc/tests/Test.cpp

rpc_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_tests_:.cpp=.o))

