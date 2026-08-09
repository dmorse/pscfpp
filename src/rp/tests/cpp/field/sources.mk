rpc_tests_field_=rpc/tests/field/Test.cpp

rpc_tests_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_tests_field_:.cpp=.o))

