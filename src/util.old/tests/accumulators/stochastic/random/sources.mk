util_tests_accumulators_stochastic_random_ = \
    util/tests/accumulators/stochastic/random/RandomTest.cc 

util_tests_accumulators_stochastic_random_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_accumulators_stochastic_random_))
util_tests_accumulators_stochastic_random_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_accumulators_stochastic_random_:.cc=.o))

