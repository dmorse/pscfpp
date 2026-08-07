# ------------------------------------------------------------------
#  Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/cpc/sources.mk
include $(SRC_DIR)/cpc/tests/sources.mk
include $(SRC_DIR)/cpc/patterns.mk

TEST=cpc/tests/Test

all: $(cpc_tests_OBJS) $(BLD_DIR)/$(TEST)

includes:
	echo $(INCLUDES)

run: $(cpc_tests_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/cpc/tests/ 

run-quiet: $(cpc_tests_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/cpc/tests/  > log
	@echo `grep failed log` ", "\
              `grep successful log` "in cpc/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	-cd $(SRC_DIR)/cpc/tests/solvers; $(MAKE) clean-outputs

clean:
	rm -f $(BLD_DIR)/$(TEST) 
	rm -f $(cpc_tests_OBJS) $(cpc_tests_OBJS:.o=.d)
	rm -f log count 
	-cd $(SRC_DIR)/cpc/tests/solvers; $(MAKE) clean
	rm -f $(SRC_DIR)/cpc/tests/*/Test
	rm -f $(SRC_DIR)/cpc/tests/*.o $(SRC_DIR)/cpc/tests/*.d
	rm -f $(SRC_DIR)/cpc/tests/*/*.o $(SRC_DIR)/cpc/tests/*/*.d

# ------------------------------------------------------------------
#  Include dependency files

-include $(cpc_tests_OBJS:.o=.d)
-include $(cpc_OBJS:.o=.d)
-include $(prdc_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
