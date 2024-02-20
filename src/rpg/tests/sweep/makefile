BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/rpg/include.mk
include $(SRC_DIR)/rpg/tests/sweep/sources.mk

TEST=rpg/tests/sweep/Test

all: $(rpg_tests_sweep_OBJS) $(BLD_DIR)/$(TEST)

includes:
	@echo $(INCLUDES)

libs:
	@echo $(LIBS)

run: $(rpg_tests_sweep_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpg/tests/ 

run-quiet: $(rpg_tests_sweep_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpg/tests/  > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rpg/tests/sweep/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	find out -type f \( ! -iname ".*" \) -exec rm -f {} \;

clean:
	rm -f $(rpg_tests_sweep_OBJS) $(rpg_tests_sweep_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	$(MAKE) clean-outputs

-include $(rpg_tests_sweep_OBJS:.o=.d)
-include $(rpg_tests_sweep_OBJS:.o=.d)
