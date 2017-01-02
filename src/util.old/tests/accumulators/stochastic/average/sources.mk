util_tests_accumulators_stochastic_average_ = \
    util/tests/accumulators/stochastic/average/AverageTest.cc 

util_tests_accumulators_stochastic_average_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_accumulators_stochastic_average_))
util_tests_accumulators_stochastic_average_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_accumulators_stochastic_average_:.cc=.o))

