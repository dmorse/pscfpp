SRC_DIR_REL =../../../../..

include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include sources.mk

all: $(util_tests_accumulators_stochastic_average_OBJS) AverageTest

clean:
	rm -f $(util_tests_accumulators_stochastic_average_OBJS) 
	rm -f $(util_tests_accumulators_stochastic_average_OBJS:.o=.d)
	rm -f AverageTest

clean-deps:
	rm -f $(util_tests_accumulators_stochastic_average_OBJS:.o=.d)

AverageTest: AverageTest.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o AverageTest AverageTest.o \
	       $(util_LIB)

-include $(util_tests_accumulators_stochastic_average_OBJS:.o=.d)

