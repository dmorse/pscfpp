BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/pscf/include.mk
include $(SRC_DIR)/pscf/tests/cpu/sources.mk

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
	rm -f $(pscf_tests_cpu_OBJS) $(pscf_tests_cpu_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST)
	$(MAKE) clean-outputs

-include $(pscf_tests_cpu_OBJS:.o=.d)
-include $(pscf_tests_cpu_OBJS:.o=.d)
