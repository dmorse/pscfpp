SRC_DIR_REL =../../../../..

include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include sources.mk

all: $(util_tests_accumulators_stochastic_symmTensorAverage_OBJS) SymmTensorAverageTest

clean:
	rm -f $(util_tests_accumulators_stochastic_symmTensorAverage_OBJS) 
	rm -f $(util_tests_accumulators_stochastic_symmTensorAverage_OBJS:.o=.d)
	rm -f SymmTensorAverageTest

clean-deps:
	rm -f $(util_tests_accumulators_stochastic_symmTensorAverage_OBJS:.o=.d)

SymmTensorAverageTest: SymmTensorAverageTest.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o SymmTensorAverageTest SymmTensorAverageTest.o \
	       $(util_LIB)

-include $(util_tests_accumulators_stochastic_symmTensorAverage_OBJS:.o=.d)

