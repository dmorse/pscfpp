r1d_tests_=r1d/tests/Test.cpp

r1d_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(r1d_tests_:.cpp=.o))

