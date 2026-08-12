# -------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/cpc/sources.mk
include $(SRC_DIR)/cpc/tests/field/sources.mk
include $(SRC_DIR)/cpc/patterns.mk

# -------------------------------------------------------------------------
# Main targets

TEST=cpc/tests/field/Test

all: $(cpc_tests_field_OBJS) $(BLD_DIR)/$(TEST)

run: $(cpc_tests_field_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/cpc/tests/ 

run-quiet: $(cpc_tests_field_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/cpc/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in cpc/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count out/*

clean:
	rm -f $(BLD_DIR)/$(TEST)
	rm -f $(cpc_tests_field_OBJS) $(cpc_tests_field_OBJS:.o=.d)
	$(MAKE) clean-outputs

.PHONY: all includes defines libs run run-quiet clean-outputs clean

# -------------------------------------------------------------------------
# Include dependency files

-include $(cpc_tests_field_OBJS:.o=.d)
-include $(cpc_OBJS:.o=.d)
-include $(prdc_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
