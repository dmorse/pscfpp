util_tests_accumulators_stochastic_tensorAverage_ = \
    util/tests/accumulators/stochastic/tensorAverage/TensorAverageTest.cc 

util_tests_accumulators_stochastic_tensorAverage_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_accumulators_stochastic_tensorAverage_))
util_tests_accumulators_stochastic_tensorAverage_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_accumulators_stochastic_tensorAverage_:.cc=.o))

