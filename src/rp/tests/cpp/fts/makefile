# -----------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/rp/sources.mk
include $(SRC_DIR)/rp/tests/cpp/fts/sources.mk
include $(SRC_DIR)/rp/patterns.mk

# -----------------------------------------------------------------------------
# Main targets

TEST=rp/tests/cpp/fts/Test

all: $(rpc_tests_fts_OBJS) $(BLD_DIR)/$(TEST)

run: $(rpc_tests_fts_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rp/tests/cpp/fts/

run-quiet: $(rpc_tests_fts_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rp/tests/cpp/fts/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in rp/tests/cpp/fts/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	rm -f out/*

clean:
	rm -f $(BLD_DIR)/$(TEST)
	rm -f $(rpc_tests_fts_OBJS) $(rpc_tests_fts_OBJS:.o=.d)
	$(MAKE) clean-outputs

# -----------------------------------------------------------------------------
# Include dependency files

-include $(rpc_tests_fts_OBJS:.o=.d)
-include $(rp_DEPS)
-include $(prdc_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
