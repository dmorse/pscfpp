BLD_DIR_REL =../../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/util/tests/accumulators/unit/sources.mk

TEST=$(BLD_DIR)/util/tests/accumulators/unit/Test

all: $(util_tests_accumulators_unit_OBJS) $(TEST)

clean:
	rm -f $(util_tests_accumulators_unit_OBJS) $(util_tests_accumulators_unit_OBJS:.o=.d)
	rm -f $(TEST) 
	$(MAKE) clean-outputs

clean-outputs:
	rm -f binary

-include $(util_tests_accumulators_unit_OBJS:.o=.d)

