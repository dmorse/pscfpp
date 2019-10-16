pspc_tests_=pspc/tests/Test.cc

pspc_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspc_tests_))
pspc_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_tests_:.cc=.o))

