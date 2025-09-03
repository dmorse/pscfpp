# ------------------------------------------------------------------
#  Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/util/config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/rpc/sources.mk
include $(SRC_DIR)/rpc/tests/sources.mk
include $(SRC_DIR)/rpc/patterns.mk

TEST=rpc/tests/Test

all: $(rpc_tests_OBJS) $(BLD_DIR)/$(TEST)

includes:
	echo $(INCLUDES)

run: $(rpc_tests_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpc/tests/ 

run-quiet: $(rpc_tests_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpc/tests/  > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rpc/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	-cd $(SRC_DIR)/rpc/tests/field; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rpc/tests/environment; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rpc/tests/solvers; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rpc/tests/system; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rpc/tests/sweep; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rpc/tests/fts; $(MAKE) clean-outputs

clean:
	rm -f $(BLD_DIR)/$(TEST) 
	rm -f $(rpc_tests_OBJS) $(rpc_tests_OBJS:.o=.d)
	rm -f log count 
	-cd $(SRC_DIR)/rpc/tests/field; $(MAKE) clean
	-cd $(SRC_DIR)/rpc/tests/environment; $(MAKE) clean
	-cd $(SRC_DIR)/rpc/tests/solvers; $(MAKE) clean
	-cd $(SRC_DIR)/rpc/tests/system; $(MAKE) clean
	-cd $(SRC_DIR)/rpc/tests/sweep; $(MAKE) clean
	-cd $(SRC_DIR)/rpc/tests/fts; $(MAKE) clean
	rm -f $(SRC_DIR)/rpc/tests/*/Test
	rm -f $(SRC_DIR)/rpc/tests/*.o $(SRC_DIR)/rpc/tests/*.d
	rm -f $(SRC_DIR)/rpc/tests/*/*.o $(SRC_DIR)/rpc/tests/*/*.d

# ------------------------------------------------------------------
#  Include dependency files

-include $(rpc_tests_OBJS:.o=.d)
-include $(rpc_OBJS:.o=.d)
-include $(prdc_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
