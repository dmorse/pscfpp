util_tests_format_=util/tests/format/Test.cc

util_tests_format_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_format_))
util_tests_format_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_format_:.cc=.o))

