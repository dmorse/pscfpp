include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/prdc/tests/crystal/sources.mk
include $(SRC_DIR)/prdc/patterns.mk

TEST=prdc/tests/crystal/Test

all: $(prdc_tests_crystal_OBJS) $(BLD_DIR)/$(TEST)

run: $(prdc_tests_crystal_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/prdc/tests/crystal/

run-quiet: $(prdc_tests_crystal_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/prdc/tests/crystal/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in prdc/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count
	rm -f out/*

clean:
	rm -f $(prdc_tests_crystal_OBJS)
	rm -f $(prdc_tests_crystal_DEPS)
	rm -f $(BLD_DIR)/$(TEST)
	$(MAKE) clean-outputs

.PHONY: all run run-quiet clean-outputs clean

-include $(util_DEPS)
-include $(pscf_DEPS)
-include $(prdc_DEPS)
-include $(prdc_tests_crystal_DEPS)
