include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk
include $(SRC_DIR)/pscf/tests/interaction/sources.mk

TEST=pscf/tests/interaction/Test

all: $(pscf_tests_interaction_OBJS) $(BLD_DIR)/$(TEST)

run: $(pscf_tests_interaction_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pscf/tests/inter/ 

run-quiet: $(pscf_tests_interaction_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/pscf/tests/inter/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pscf/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 

clean:
	rm -f $(pscf_tests_interaction_OBJS) 
	rm -f $(pscf_tests_interaction_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST)
	$(MAKE) clean-outputs

-include $(pscf_tests_interaction_OBJS:.o=.d)
-include $(pscf_DEPS)
-include $(util_DEPS)
