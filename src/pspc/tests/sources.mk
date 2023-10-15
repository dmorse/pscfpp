pspc_tests_=pspc/tests/Test.cpp

pspc_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_tests_:.cpp=.o))

