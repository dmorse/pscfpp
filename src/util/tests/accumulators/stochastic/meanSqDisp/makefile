SRC_DIR_REL =../../../../..

include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include sources.mk

all: $(util_tests_accumulators_stochastic_meanSqDisp_OBJS) MeanSqDispTest

clean:
	rm -f $(util_tests_accumulators_stochastic_meanSqDisp_OBJS) 
	rm -f $(util_tests_accumulators_stochastic_meanSqDisp_OBJS:.o=.d)
	rm -f MeanSqDispTest

clean-deps:
	rm -f $(util_tests_accumulators_stochastic_meanSqDisp_OBJS:.o=.d)

MeanSqDispTest: MeanSqDispTest.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o MeanSqDispTest MeanSqDispTest.o \
	       $(util_LIB)

-include $(util_tests_accumulators_stochastic_meanSqDisp_OBJS:.o=.d)

