include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/tests/cpu/sources.mk
include $(SRC_DIR)/prdc/patterns.mk

TEST=prdc/tests/cpu/Test

all: $(prdc_tests_cpu_OBJS) $(BLD_DIR)/$(TEST)

run: $(prdc_tests_cpu_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/prdc/tests/

run-quiet: $(prdc_tests_cpu_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/prdc/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in prdc/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count out/*

clean:
	rm -f $(prdc_tests_cpu_OBJS)
	rm -f $(prdc_tests_cpu_DEPS)
	rm -f $(BLD_DIR)/$(TEST)
	$(MAKE) clean-outputs

.PHONY: all run run-quiet clean-outputs clean

-include $(prdc_tests_cpu_DEPS)
-include $(prdc_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
