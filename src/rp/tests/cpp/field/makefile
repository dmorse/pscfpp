# -------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/rp/sources.mk
include $(SRC_DIR)/rp/patterns.mk
include $(SRC_DIR)/rp/tests/cpp/field/sources.mk

# -------------------------------------------------------------------------
# Main targets

TEST=rp/tests/cpp/field/Test

all: $(rpc_tests_field_OBJS) $(BLD_DIR)/$(TEST)

run: $(rpc_tests_field_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rp/tests/cpp/ 

run-quiet: $(rpc_tests_field_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rp/tests/cpp/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rp/tests/cpp/log" > count
	@cat count

clean-outputs:
	rm -f log count out/*

clean:
	rm -f $(BLD_DIR)/$(TEST)
	rm -f $(rpc_tests_field_OBJS) $(rpc_tests_field_OBJS:.o=.d)
	$(MAKE) clean-outputs

.PHONY: all includes defines libs run run-quiet clean-outputs clean

# -------------------------------------------------------------------------
# Include dependency files

-include $(rpc_tests_field_OBJS:.o=.d)
-include $(rp_DEPS)
-include $(prdc_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
