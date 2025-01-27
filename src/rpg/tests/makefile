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
ifeq ($(SRC_DIR),$(BLD_DIR))
	-cd system; $(MAKE) clean-outputs
	-cd field; $(MAKE) clean-outputs
	-cd fts; $(MAKE) clean-outputs
	-cd solvers; $(MAKE) clean-outputs
	-cd sweep; $(MAKE) clean-outputs
endif

clean:
	rm -f $(rpg_tests_OBJS) $(rpg_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	rm -f log count 
ifeq ($(SRC_DIR),$(BLD_DIR))
	-cd system; $(MAKE) clean
	-cd field; $(MAKE) clean
	-cd fts; $(MAKE) clean
	-cd solvers; $(MAKE) clean
	-cd sweep; $(MAKE) clean
	rm -f */Test */*.o */*.d 
endif

-include $(rpg_tests_OBJS:.o=.d)
-include $(rpg_OBJS:.o=.d)
