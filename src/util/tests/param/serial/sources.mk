util_tests_param_serial_=util/tests/param/serial/Test.cc

util_tests_param_serial_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_param_serial_))
util_tests_param_serial_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_param_serial_:.cc=.o))

