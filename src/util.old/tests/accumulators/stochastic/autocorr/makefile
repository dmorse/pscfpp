SRC_DIR_REL =../../../../..

include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk


include sources.mk

all: $(util_tests_accumulators_stochastic_autocorr_OBJS) AutoCorrTest

clean:
	rm -f $(util_tests_accumulators_stochastic_autocorr_OBJS)
	rm -f $(util_tests_accumulators_stochastic_autocorr_OBJS:.o=.d)
	rm -f AutoCorrTest 

clean-deps:
	rm -f $(util_tests_accumulators_stochastic_autocorr_OBJS:.o=.d)

AutoCorrTest: AutoCorrTest.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o AutoCorrTest AutoCorrTest.o \
	       $(util_LIB)

-include $(util_tests_accumulators_stochastic_autocorr_OBJS:.o=.d)

