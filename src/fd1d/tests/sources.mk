fd1d_tests_=fd1d/tests/Test.cc

fd1d_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(fd1d_tests_))
fd1d_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(fd1d_tests_:.cc=.o))

