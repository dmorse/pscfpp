# -------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/cpc/sources.mk
include $(SRC_DIR)/cpc/tests/solvers/sources.mk
include $(SRC_DIR)/cpc/patterns.mk

# -------------------------------------------------------------------------
# Main targets

TEST=cpc/tests/solvers/Test

all: $(cpc_tests_solvers_OBJS) $(BLD_DIR)/$(TEST)

run: $(cpc_tests_solvers_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/cpc/tests/ 

run-quiet: $(cpc_tests_solvers_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/cpc/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in cpc/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count

clean:
	rm -f $(BLD_DIR)/$(TEST)
	rm -f $(cpc_tests_solvers_OBJS) $(cpc_tests_solvers_OBJS:.o=.d)
	$(MAKE) clean-outputs

# -------------------------------------------------------------------------
# Main targets

-include $(cpc_tests_solvers_OBJS:.o=.d)
-include $(cpc_OBJS:.o=.d)
-include $(prdc_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
