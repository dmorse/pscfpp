util_tests_accumulators_stochastic_autoCorrelation_ = \
    util/tests/accumulators/stochastic/autoCorrelation/AutoCorrelationTest.cc 

util_tests_accumulators_stochastic_autoCorrelation_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_accumulators_stochastic_autoCorrelation_))
util_tests_accumulators_stochastic_autoCorrelation_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_accumulators_stochastic_autoCorrelation_:.cc=.o))

