include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/rpc/sources.mk
include $(SRC_DIR)/rpc/tests/environment/sources.mk
include $(SRC_DIR)/rpc/patterns.mk

TEST=rpc/tests/environment/Test

.PHONY: all run run-quiet clean-outputs clean

all: $(rpc_tests_environment_OBJS) $(BLD_DIR)/$(TEST)

run: $(rpc_tests_environment_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpc/tests/ 

run-quiet: $(rpc_tests_environment_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpc/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rpc/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count out/*.*
	rm -f out/sweep/*

clean:
	rm -f $(BLD_DIR)/$(TEST)
	rm -f $(rpc_tests_environment_OBJS) $(rpc_tests_environment_OBJS:.o=.d)
	$(MAKE) clean-outputs

-include $(TEST).d
-include $(rpc_DEPS)
-include $(prdc_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
