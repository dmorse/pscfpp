BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/rpc/include.mk
include $(SRC_DIR)/rpc/tests/system/sources.mk

TEST=rpc/tests/system/Test

all: $(rpc_tests_system_OBJS) $(BLD_DIR)/$(TEST)

includes:
	@echo $(INCLUDES)

libs:
	@echo $(LIBS)

run: $(rpc_tests_system_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpc/tests/system/

run-quiet: $(rpc_tests_system_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpc/tests/system/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in rpc/tests/system/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	rm -f out/*

clean:
	rm -f $(BLD_DIR)/$(TEST)
	rm -f $(rpc_tests_system_OBJS) $(rpc_tests_system_OBJS:.o=.d)
	$(MAKE) clean-outputs

-include $(rpc_tests_system_OBJS:.o=.d)
-include $(rpc_tests_system_OBJS:.o=.d)
