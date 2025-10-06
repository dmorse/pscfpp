include ../../../config.mk
include $(SRC_DIR)/rpg/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/rpg/sources.mk
include $(SRC_DIR)/rpg/tests/environment/sources.mk

TEST=rpg/tests/environment/Test

all: $(rpg_tests_environment_OBJS) $(BLD_DIR)/$(TEST)

run: $(rpg_tests_environment_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpg/tests/ 

run-quiet: $(rpg_tests_environment_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpg/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rpg/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count out/*.*
	rm -f out/sweep/*

clean:
	rm -f $(BLD_DIR)/$(TEST)
	rm -f $(rpg_tests_environment_OBJS) $(rpg_tests_environment_OBJS:.o=.d)
	$(MAKE) clean-outputs

.PHONY: all run run-quiet clean-outputs clean

-include $(util_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(prdc_OBJS:.o=.d)
-include $(rpg_OBJS:.o=.d)
-include $(rpg_tests_environment_OBJS:.o=.d)
