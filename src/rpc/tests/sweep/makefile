BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/rpc/include.mk
include $(SRC_DIR)/rpc/tests/sweep/sources.mk

TEST=rpc/tests/sweep/Test

all: $(rpc_tests_sweep_OBJS) $(BLD_DIR)/$(TEST)

includes:
	@echo $(INCLUDES)

libs:
	@echo $(LIBS)

run: $(rpc_tests_sweep_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpc/tests/ 

run-quiet: $(rpc_tests_sweep_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpc/tests/  > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rpc/tests/sweep/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	rm -f out/*/*
	rm -f out/test*

#find out -type f \( ! -iname ".*" \) -exec rm -f {} \;

clean:
	rm -f $(BLD_DIR)/$(TEST)
	rm -f $(rpc_tests_sweep_OBJS) $(rpc_tests_sweep_OBJS:.o=.d)
	$(MAKE) clean-outputs

-include $(rpc_tests_sweep_OBJS:.o=.d)
-include $(rpc_tests_sweep_OBJS:.o=.d)
