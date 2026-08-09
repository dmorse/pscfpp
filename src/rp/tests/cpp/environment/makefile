include ../../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/rp/sources.mk
include $(SRC_DIR)/rp/tests/cpp/environment/sources.mk
include $(SRC_DIR)/rp/patterns.mk

TEST=rp/tests/cpp/environment/Test

.PHONY: all run run-quiet clean-outputs clean

all: $(rpc_tests_environment_OBJS) $(BLD_DIR)/$(TEST)

run: $(rpc_tests_environment_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rp/tests/cpp/

run-quiet: $(rpc_tests_environment_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rp/tests/cpp/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rp/tests/cpp/log" > count
	@cat count

clean-outputs:
	rm -f log count out/*.*
	rm -f out/sweep/*

clean:
	rm -f $(BLD_DIR)/$(TEST)
	rm -f $(rpc_tests_environment_OBJS) $(rpc_tests_environment_OBJS:.o=.d)
	$(MAKE) clean-outputs

-include $(TEST).d
-include $(rp_DEPS)
-include $(prdc_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
