BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/rpc/include.mk
include $(SRC_DIR)/rpc/tests/simulate/sources.mk

TEST=rpc/tests/simulate/Test

all: $(rpc_tests_simulate_OBJS) $(BLD_DIR)/$(TEST)

includes:
	@echo $(INCLUDES)

libs:
	@echo $(LIBS)

run: $(rpc_tests_simulate_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpc/tests/simulate/

run-quiet: $(rpc_tests_simulate_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpc/tests/simulate/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in rpc/tests/simulate/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	rm -f out/*

clean:
	rm -f $(rpc_tests_simulate_OBJS) $(rpc_tests_simulate_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	$(MAKE) clean-outputs

-include $(rpc_tests_simulate_OBJS:.o=.d)
-include $(rpc_tests_simulate_OBJS:.o=.d)
