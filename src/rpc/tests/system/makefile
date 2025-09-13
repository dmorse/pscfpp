# -----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/rpc/sources.mk
include $(SRC_DIR)/rpc/tests/system/sources.mk
include $(SRC_DIR)/rpc/patterns.mk

# -----------------------------------------------------------------------
# Main targets

TEST=rpc/tests/system/Test

all: $(rpc_tests_system_OBJS) $(BLD_DIR)/$(TEST)

run: $(rpc_tests_system_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpc/tests/system/

run-quiet: $(rpc_tests_system_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpc/tests/system/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in rpc/tests/system/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	rm -f out/*

clean:
	rm -f $(BLD_DIR)/$(TEST)
	rm -f $(rpc_tests_system_OBJS) $(rpc_tests_system_OBJS:.o=.d)
	$(MAKE) clean-outputs

# -----------------------------------------------------------------------
# Include dependency files

-include $(rpc_tests_system_OBJS:.o=.d)
-include $(rpc_OBJS:.o=.d)
-include $(prdc_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
