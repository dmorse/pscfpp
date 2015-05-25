util_tests_accumulators_stochastic_averageStage_ = \
    util/tests/accumulators/stochastic/averageStage/AverageStageTest.cc 

util_tests_accumulators_stochastic_averageStage_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_accumulators_stochastic_averageStage_))
util_tests_accumulators_stochastic_averageStage_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_accumulators_stochastic_averageStage_:.cc=.o))

