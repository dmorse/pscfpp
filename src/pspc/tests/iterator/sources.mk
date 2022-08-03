pspc_tests_iterator_=pspc/tests/iterator/Test.cc

pspc_tests_iterator_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspc_tests_iterator_))
pspc_tests_iterator_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_tests_iterator_:.cc=.o))
