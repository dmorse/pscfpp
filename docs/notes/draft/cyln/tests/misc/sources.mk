cyln_tests_misc_=cyln/tests/misc/Test.cc

cyln_tests_misc_SRCS=\
     $(addprefix $(SRC_DIR)/, $(cyln_tests_misc_))
cyln_tests_misc_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cyln_tests_misc_:.cc=.o))

