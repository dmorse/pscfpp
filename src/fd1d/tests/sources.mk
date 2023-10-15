fd1d_tests_=fd1d/tests/Test.cpp

fd1d_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(fd1d_tests_:.cpp=.o))

