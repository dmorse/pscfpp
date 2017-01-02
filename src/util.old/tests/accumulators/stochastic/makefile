SRC_DIR_REL =../../../..

include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk

include sources.mk

all: $(tests_util_accumulators_stochastic_OBJS)
	cd autocorr; $(MAKE) all
	cd autoCorrelation; $(MAKE) all
	cd average;  $(MAKE) all
	cd averageStage; $(MAKE) all
	cd meanSqDisp; $(MAKE) all
	cd random; $(MAKE) all

clean:
	cd autocorr; $(MAKE) clean
	cd autoCorrelation; $(MAKE) clean
	cd average;  $(MAKE) clean
	cd averageStage; $(MAKE) clean
	cd meanSqDisp; $(MAKE) clean
	cd random; $(MAKE) clean
	rm -f $(tests_util_accumulators_stochastic_OBJS) $(tests_util_accumulators_stochastic_OBJS:.o=.d)

clean-deps:
	rm -f $(tests_util_accumulators_stochastic_OBJS:.o=.d)

-include $(tests_util_accumulators_stochastic_OBJS:.o=.d)

