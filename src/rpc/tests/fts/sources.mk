rpc_tests_fts_=rpc/tests/fts/Test.cpp

rpc_tests_fts_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_tests_fts_:.cpp=.o))

