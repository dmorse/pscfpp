# -------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/rp/sources.mk
include $(SRC_DIR)/rp/tests/cpp/sweep/sources.mk
include $(SRC_DIR)/rp/patterns.mk

# -------------------------------------------------------------------------
# Main targets

TEST=rp/tests/cpp/sweep/Test

all: $(rpc_tests_sweep_OBJS) $(BLD_DIR)/$(TEST)

run: $(rpc_tests_sweep_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rp/tests/cpp/ 

run-quiet: $(rpc_tests_sweep_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/rp/tests/cpp/  > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rp/tests/cpp/sweep/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	rm -f out/*/*
	rm -f out/test*

#find out -type f \( ! -iname ".*" \) -exec rm -f {} \;

clean:
	rm -f $(BLD_DIR)/$(TEST)
	rm -f $(rpc_tests_sweep_OBJS) $(rpc_tests_sweep_OBJS:.o=.d)
	$(MAKE) clean-outputs

# -------------------------------------------------------------------------
# Include dependency files

-include $(rpc_tests_sweep_OBJS:.o=.d)
-include $(rp_DEPS)
-include $(prdc_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
