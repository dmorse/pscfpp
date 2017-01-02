util_tests_misc_=util/tests/misc/Test.cc

util_tests_misc_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_misc_))
util_tests_misc_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_misc_:.cc=.o))

