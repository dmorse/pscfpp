BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/pscf/include.mk
include $(SRC_DIR)/pscf/tests/homogeneous/sources.mk

TEST=pscf/tests/homogeneous/Test

all: $(pscf_tests_homogeneous_OBJS) $(BLD_DIR)/$(TEST)

includes:
	echo $(INCLUDES)

run: $(pscf_tests_homogeneous_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pscf/tests/homogeneous/ 

run-quiet: $(pscf_tests_homogeneous_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pscf/tests/homogeneous/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pscf/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 

clean:
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	rm -f $(pscf_tests_homogeneous_OBJS) $(pscf_tests_homogeneous_OBJS:.o=.d)
	$(MAKE) clean-outputs

-include $(pscf_tests_homogeneous_OBJS:.o=.d)
-include $(pscf_tests_homogeneous_OBJS:.o=.d)
