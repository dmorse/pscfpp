cpc_tests_=cpc/tests/Test.cpp

cpc_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cpc_tests_:.cpp=.o))

