
include ../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/r1d/sources.mk
include $(SRC_DIR)/r1d/tests/sources.mk
include $(SRC_DIR)/r1d/patterns.mk

TEST=r1d/tests/Test

all: $(r1d_tests_OBJS) $(BLD_DIR)/$(TEST)

run: $(r1d_tests_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/r1d/tests/ 

run-quiet: $(r1d_tests_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/r1d/tests/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in r1d/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	rm -f $(SRC_DIR)/r1d/tests/out/*

clean:
	rm -f $(r1d_tests_OBJS) $(r1d_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	$(MAKE) clean-outputs

-include $(r1d_tests_OBJS:.o=.d)
-include $(r1d_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
