pssp_tests_iterator_=pssp/tests/iterator/Test.cc

pssp_tests_iterator_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_tests_iterator_))
pssp_tests_iterator_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_tests_iterator_:.cc=.o))

