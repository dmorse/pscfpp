SRC_DIR_REL =../../../../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include sources.mk

all: $(util_tests_accumulators_stochastic_averageStage_OBJS) AverageStageTest

clean:
	rm -f $(util_tests_accumulators_stochastic_averageStage_OBJS) 
	rm -f $(util_tests_accumulators_stochastic_averageStage_OBJS:.o=.d)
	rm -f AverageStageTest

clean-deps:
	rm -f $(util_tests_accumulators_stochastic_averageStage_OBJS:.o=.d)

AverageStageTest: AverageStageTest.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o AverageStageTest AverageStageTest.o \
	       $(util_LIB)

-include $(util_tests_accumulators_stochastic_averageStage_OBJS:.o=.d)

