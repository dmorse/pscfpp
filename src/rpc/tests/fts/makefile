# -----------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/util/config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/rpc/sources.mk
include $(SRC_DIR)/rpc/tests/fts/sources.mk
include $(SRC_DIR)/rpc/patterns.mk

# -----------------------------------------------------------------------------
# Main targets

TEST=rpc/tests/fts/Test

all: $(rpc_tests_fts_OBJS) $(BLD_DIR)/$(TEST)

run: $(rpc_tests_fts_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpc/tests/fts/

run-quiet: $(rpc_tests_fts_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpc/tests/fts/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in rpc/tests/fts/log" > count
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
-include $(rpc_OBJS:.o=.d)
-include $(prdc_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
