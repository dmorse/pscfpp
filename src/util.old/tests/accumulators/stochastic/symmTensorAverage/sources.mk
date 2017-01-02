util_tests_accumulators_stochastic_symmTensorAverage_ = \
    util/tests/accumulators/stochastic/symmTensorAverage/SymmTensorAverageTest.cc 

util_tests_accumulators_stochastic_symmTensorAverage_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_accumulators_stochastic_symmTensorAverage_))
util_tests_accumulators_stochastic_symmTensorAverage_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_accumulators_stochastic_symmTensorAverage_:.cc=.o))

