rpc_tests_system_=rpc/tests/system/Test.cpp

rpc_tests_system_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_tests_system_:.cpp=.o))

