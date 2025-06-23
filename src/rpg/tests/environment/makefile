BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/rpg/include.mk
include $(SRC_DIR)/rpg/tests/environment/sources.mk

TEST=rpg/tests/environment/Test

.PHONY: all includes defines libs run run-quiet clean-outputs clean

all: $(rpg_tests_environment_OBJS) $(BLD_DIR)/$(TEST)

includes:
	@echo $(INCLUDES)

defines:
	@echo $(DEFINES)

libs:
	@echo $(LIBS)

run: $(rpg_tests_environment_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpg/tests/ 

run-quiet: $(rpg_tests_environment_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpg/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rpg/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count out/*.*
	rm -f out/sweep/*

clean:
	rm -f $(BLD_DIR)/$(TEST)
	rm -f $(rpg_tests_environment_OBJS) $(rpg_tests_environment_OBJS:.o=.d)
	$(MAKE) clean-outputs

-include $(rpg_OBJS:.o=.d)
-include $(rpg_tests_environment_OBJS:.o=.d)
