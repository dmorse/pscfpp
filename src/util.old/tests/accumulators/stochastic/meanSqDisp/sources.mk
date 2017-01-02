util_tests_accumulators_stochastic_meanSqDisp_ = \
    util/tests/accumulators/stochastic/meanSqDisp/MeanSqDispTest.cc 

util_tests_accumulators_stochastic_meanSqDisp_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_accumulators_stochastic_meanSqDisp_))
util_tests_accumulators_stochastic_meanSqDisp_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_accumulators_stochastic_meanSqDisp_:.cc=.o))

