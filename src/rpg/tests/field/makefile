include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/rpg/sources.mk
include $(SRC_DIR)/rpg/tests/field/sources.mk
include $(SRC_DIR)/rpg/patterns.mk

TEST=rpg/tests/field/Test

all: $(rpg_tests_field_OBJS) $(BLD_DIR)/$(TEST)

run: $(rpg_tests_field_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpg/tests/field/

run-quiet: $(rpg_tests_field_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpg/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rpg/tests/field/log" > count
	@cat count

clean-outputs:
	rm -f log count out/*

clean:
	rm -f $(rpg_tests_field_OBJS) $(rpg_tests_field_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	$(MAKE) clean-outputs

.PHONY: all run run-quiet clean-outputs clean

-include $(util_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(prdc_OBJS:.o=.d)
-include $(rpg_OBJS:.o=.d)
-include $(rpg_tests_field_OBJS:.o=.d)
