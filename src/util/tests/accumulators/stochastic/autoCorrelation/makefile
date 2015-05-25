SRC_DIR_REL =../../../../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include sources.mk

all: $(util_tests_accumulators_stochastic_autoCorrelation_OBJS) AutoCorrelationTest

clean:
	rm -f $(util_tests_accumulators_stochastic_autoCorrelation_OBJS) 
	rm -f $(util_tests_accumulators_stochastic_autoCorrelation_OBJS:.o=.d)
	rm -f AutoCorrelationTest
	rm -f binary

clean-deps:
	rm -f $(util_tests_accumulators_stochastic_autoCorrelation_OBJS:.o=.d)

AutoCorrelationTest: AutoCorrelationTest.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o AutoCorrelationTest AutoCorrelationTest.o \
	       $(util_LIB)

-include $(util_tests_accumulators_stochastic_autoCorrelation_OBJS:.o=.d)

