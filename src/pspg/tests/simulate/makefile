BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/pspg/include.mk
include $(SRC_DIR)/pspg/tests/simulate/sources.mk

TEST=pspg/tests/simulate/Test

all: $(pspg_tests_simulate_OBJS) $(BLD_DIR)/$(TEST)

includes:
	@echo $(INCLUDES)

libs:
	@echo $(LIBS)

run: $(pspg_tests_simulate_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pspg/tests/simulate/

run-quiet: $(pspg_tests_simulate_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pspg/tests/simulate/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in pspg/tests/simulate/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	rm -f out/*

clean:
	rm -f $(pspg_tests_simulate_OBJS) $(pspg_tests_simulate_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	$(MAKE) clean-outputs

-include $(pspg_tests_simulate_OBJS:.o=.d)
-include $(pspg_tests_simulate_OBJS:.o=.d)
