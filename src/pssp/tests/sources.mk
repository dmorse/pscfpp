pssp_tests_=pssp/tests/Test.cc

pssp_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_tests_))
pssp_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_tests_:.cc=.o))

