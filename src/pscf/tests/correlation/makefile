include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(BLD_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/tests/correlation/sources.mk
include $(BLD_DIR)/pscf/patterns.mk

TEST=pscf/tests/correlation/Test

all: $(pscf_tests_correlation_OBJS) $(BLD_DIR)/$(TEST)

includes:
	echo $(INCLUDES)

run: $(pscf_tests_correlation_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pscf/tests/correlation/ 

run-quiet: $(pscf_tests_correlation_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pscf/tests/correlation/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pscf/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 

clean:
	rm -f $(pscf_tests_correlation_OBJS)
	rm -f $(pscf_tests_correlation_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST)
	$(MAKE) clean-outputs

-include $(pscf_tests_correlation_OBJS:.o=.d)
-include $(pscf_DEPS)
-include $(util_DEPS)
