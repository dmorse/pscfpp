util_tests_boundary_=util/tests/boundary/Test.cc

util_tests_boundary_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_boundary_))
util_tests_boundary_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_boundary_:.cc=.o))

