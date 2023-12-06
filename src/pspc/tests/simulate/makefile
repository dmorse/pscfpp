BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/pspc/include.mk
include $(SRC_DIR)/pspc/tests/simulate/sources.mk

TEST=pspc/tests/simulate/Test

all: $(pspc_tests_simulate_OBJS) $(BLD_DIR)/$(TEST)

includes:
	@echo $(INCLUDES)

libs:
	@echo $(LIBS)

run: $(pspc_tests_simulate_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pspc/tests/simulate/

run-quiet: $(pspc_tests_simulate_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pspc/tests/simulate/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in pspc/tests/simulate/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	rm -f out/*

clean:
	rm -f $(pspc_tests_simulate_OBJS) $(pspc_tests_simulate_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	$(MAKE) clean-outputs

-include $(pspc_tests_simulate_OBJS:.o=.d)
-include $(pspc_tests_simulate_OBJS:.o=.d)
