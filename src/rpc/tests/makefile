BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/rpc/include.mk
include $(SRC_DIR)/rpc/tests/sources.mk

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
	-cd $(SRC_DIR)/rpc/tests/iterator; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rpc/tests/solvers; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rpc/tests/system; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rpc/tests/sweep; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rpc/tests/fts; $(MAKE) clean-outputs

clean:
	rm -f $(BLD_DIR)/$(TEST) 
	rm -f $(rpc_tests_OBJS) $(rpc_tests_OBJS:.o=.d)
	rm -f log count 
	-cd $(SRC_DIR)/rpc/tests/field; $(MAKE) clean
	-cd $(SRC_DIR)/rpc/tests/iterator; $(MAKE) clean
	-cd $(SRC_DIR)/rpc/tests/solvers; $(MAKE) clean
	-cd $(SRC_DIR)/rpc/tests/system; $(MAKE) clean
	-cd $(SRC_DIR)/rpc/tests/sweep; $(MAKE) clean
	-cd $(SRC_DIR)/rpc/tests/fts; $(MAKE) clean
	rm -f $(SRC_DIR)/rpc/tests/*/Test
	rm -f $(SRC_DIR)/rpc/tests/*.o $(SRC_DIR)/rpc/tests/*.d
	rm -f $(SRC_DIR)/rpc/tests/*/*.o $(SRC_DIR)/rpc/tests/*/*.d

-include $(rpc_tests_OBJS:.o=.d)
-include $(rpc_OBJS:.o=.d)
