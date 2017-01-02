util_tests_accumulators_unit_= \
    util/tests/accumulators/unit/Test.cc

util_tests_accumulators_unit_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_accumulators_unit_))
util_tests_accumulators_unit_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_accumulators_unit_:.cc=.o))

