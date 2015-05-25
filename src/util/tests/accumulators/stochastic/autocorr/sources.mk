util_tests_accumulators_stochastic_autocorr_ = \
    util/tests/accumulators/stochastic/autocorr/AutoCorrTest.cc 

util_tests_accumulators_stochastic_autocorr_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_accumulators_stochastic_autocorr_))
util_tests_accumulators_stochastic_autocorr_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_accumulators_stochastic_autocorr_:.cc=.o))

