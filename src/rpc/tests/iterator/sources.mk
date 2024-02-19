rpc_tests_iterator_=rpc/tests/iterator/Test.cpp

rpc_tests_iterator_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_tests_iterator_:.cpp=.o))
