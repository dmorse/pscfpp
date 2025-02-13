BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/rpg/include.mk
include $(SRC_DIR)/rpg/tests/sources.mk

TEST=rpg/tests/Test

all: $(rpg_tests_OBJS) $(BLD_DIR)/$(TEST)

includes:
	echo $(INCLUDES)

run: $(rpg_tests_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpg/tests/ 

run-quiet: $(rpg_tests_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpg/tests/  > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rpg/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	-cd $(SRC_DIR)/rpg/tests/field; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rpg/tests/iterator; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rpg/tests/solvers; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rpg/tests/system; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rpg/tests/sweep; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/rpg/tests/fts; $(MAKE) clean-outputs

clean:
	rm -f $(rpg_tests_OBJS) $(rpg_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	rm -f log count 
	-cd $(SRC_DIR)/rpg/tests/field; $(MAKE) clean
	-cd $(SRC_DIR)/rpg/tests/iterator; $(MAKE) clean
	-cd $(SRC_DIR)/rpg/tests/solvers; $(MAKE) clean
	-cd $(SRC_DIR)/rpg/tests/system; $(MAKE) clean
	-cd $(SRC_DIR)/rpg/tests/sweep; $(MAKE) clean
	-cd $(SRC_DIR)/rpg/tests/fts; $(MAKE) clean
	rm -f $(SRC_DIR)/rpg/tests/*/Test
	rm -f $(SRC_DIR)/rpg/tests/*.o $(SRC_DIR)/rpg/tests/*.d
	rm -f $(SRC_DIR)/rpg/tests/*/*.o $(SRC_DIR)/rpg/tests/*/*.d

-include $(rpg_tests_OBJS:.o=.d)
-include $(rpg_OBJS:.o=.d)
