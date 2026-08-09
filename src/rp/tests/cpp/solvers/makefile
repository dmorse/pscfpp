# -------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/rp/sources.mk
include $(SRC_DIR)/rp/tests/cpp/solvers/sources.mk
include $(SRC_DIR)/rp/patterns.mk

# -------------------------------------------------------------------------
# Main targets

TEST=rp/tests/cpp/solvers/Test

all: $(rpc_tests_solvers_OBJS) $(BLD_DIR)/$(TEST)

run: $(rpc_tests_solvers_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rp/tests/cpp/ 

run-quiet: $(rpc_tests_solvers_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rp/tests/cpp/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rp/tests/cpp/log" > count
	@cat count

clean-outputs:
	rm -f log count

clean:
	rm -f $(BLD_DIR)/$(TEST)
	rm -f $(rpc_tests_solvers_OBJS) $(rpc_tests_solvers_OBJS:.o=.d)
	$(MAKE) clean-outputs

# -------------------------------------------------------------------------
# Main targets

-include $(rpc_tests_solvers_OBJS:.o=.d)
-include $(rp_DEPS)
-include $(prdc_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
