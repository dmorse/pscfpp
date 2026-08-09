# ------------------------------------------------------------------
#  Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/rp/sources.mk
include $(SRC_DIR)/rp/patterns.mk
include $(SRC_DIR)/rp/tests/cpp/sources.mk

TEST=rp/tests/cpp/Test

all: $(rpc_tests_OBJS) $(BLD_DIR)/$(TEST)

includes:
	echo $(INCLUDES)

run: $(rpc_tests_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rp/tests/cpp/

run-quiet: $(rpc_tests_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/rp/tests/cpp/  > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rp/tests/cpp/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	-cd $(SRC_DIR)/rp/tests/cpp/field; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rp/tests/cpp/environment; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rp/tests/cpp/solvers; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rp/tests/cpp/system; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rp/tests/cpp/sweep; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rp/tests/cpp/fts; $(MAKE) clean-outputs

clean:
	rm -f $(BLD_DIR)/$(TEST) 
	rm -f $(rpc_tests_OBJS) $(rpc_tests_OBJS:.o=.d)
	rm -f log count 
	-cd $(SRC_DIR)/rp/tests/cpp/field; $(MAKE) clean
	-cd $(SRC_DIR)/rp/tests/cpp/environment; $(MAKE) clean
	-cd $(SRC_DIR)/rp/tests/cpp/solvers; $(MAKE) clean
	-cd $(SRC_DIR)/rp/tests/cpp/system; $(MAKE) clean
	-cd $(SRC_DIR)/rp/tests/cpp/sweep; $(MAKE) clean
	-cd $(SRC_DIR)/rp/tests/cpp/fts; $(MAKE) clean
	rm -f $(SRC_DIR)/rp/tests/cpp/*/Test
	rm -f $(SRC_DIR)/rp/tests/cpp/*.o $(SRC_DIR)/rp/tests/cpp/*.d
	rm -f $(SRC_DIR)/rp/tests/cpp/*/*.o $(SRC_DIR)/rp/tests/cpp/*/*.d

# ------------------------------------------------------------------
#  Include dependency files

-include $(rpc_tests_OBJS:.o=.d)
-include $(rpc_DEPS)
-include $(prdc_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
