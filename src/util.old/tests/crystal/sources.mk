util_tests_crystal_=util/tests/crystal/Test.cc

util_tests_crystal_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_crystal_))
util_tests_crystal_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_crystal_:.cc=.o))

