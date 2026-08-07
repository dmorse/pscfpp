include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/tests/cpu/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

TEST=pscf/tests/cpu/Test

all: $(pscf_tests_cpu_OBJS) $(BLD_DIR)/$(TEST)

includes:
	echo $(INCLUDES)

run: $(pscf_tests_cpu_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pscf/tests/cpu 

run-quiet: $(pscf_tests_cpu_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/pscf/tests/cpu > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pscf/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 

clean:
	rm -f $(pscf_tests_cpu_OBJS)
	rm -f $(pscf_tests_cpu_DEPS)
	rm -f $(BLD_DIR)/$(TEST)
	$(MAKE) clean-outputs

-include $(pscf_tests_cpu_DEPS)
-include $(pscf_util_DEPS)
-include $(pscf_pscf_DEPS)
