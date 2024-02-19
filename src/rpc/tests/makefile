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
ifeq ($(SRC_DIR),$(BLD_DIR))
	-cd field; $(MAKE) clean-outputs
	-cd iterator; $(MAKE) clean-outputs
	-cd solvers; $(MAKE) clean-outputs
	-cd system; $(MAKE) clean-outputs
	-cd sweep; $(MAKE) clean-outputs
else
	-cd $(SRC_DIR)/rpc/tests/field; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rpc/tests/iterator; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rpc/tests/solvers; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rpc/tests/system; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rpc/tests/sweep; $(MAKE) clean-outputs
endif

clean:
	rm -f $(rpc_tests_OBJS) $(rpc_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	rm -f log count 
ifeq ($(SRC_DIR),$(BLD_DIR))
	-cd field; $(MAKE) clean
	-cd iterator; $(MAKE) clean
	-cd solvers; $(MAKE) clean
	-cd system; $(MAKE) clean
	-cd sweep; $(MAKE) clean
	rm -f */Test */*.o */*.d
else
	-cd $(SRC_DIR)/rpc/tests/field; $(MAKE) clean
	-cd $(SRC_DIR)/rpc/tests/iterator; $(MAKE) clean
	-cd $(SRC_DIR)/rpc/tests/solvers; $(MAKE) clean
	-cd $(SRC_DIR)/rpc/tests/system; $(MAKE) clean
	-cd $(SRC_DIR)/rpc/tests/sweep; $(MAKE) clean
endif

-include $(rpc_tests_OBJS:.o=.d)
-include $(rpc_OBJS:.o=.d)
