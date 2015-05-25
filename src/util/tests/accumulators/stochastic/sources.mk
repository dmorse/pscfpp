include $(SRC_DIR)/util/tests/accumulators/stochastic/autocorr/sources.mk
include $(SRC_DIR)/util/tests/accumulators/stochastic/autoCorrelation/sources.mk
include $(SRC_DIR)/util/tests/accumulators/stochastic/average/sources.mk
include $(SRC_DIR)/util/tests/accumulators/stochastic/averageStage/sources.mk
include $(SRC_DIR)/util/tests/accumulators/stochastic/meanSqDisp/sources.mk
include $(SRC_DIR)/util/tests/accumulators/stochastic/random/sources.mk

util_tests_accumulators_stochastic_=\
    $(util_tests_accumulators_stochastic_autocorr_) \
    $(util_tests_accumulators_stochastic_autoCorrelation_) \
    $(util_tests_accumulators_stochastic_average_) \
    $(util_tests_accumulators_stochastic_averageStage_) \
    $(util_tests_accumulators_stochastic_meanSqDisp_) \
    $(util_tests_accumulators_stochastic_random_) 

util_tests_accumulators_stochastic_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_accumulators_stochastic_))
util_tests_accumulators_stochastic_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_accumulators_stochastic_:.cc=.o))

