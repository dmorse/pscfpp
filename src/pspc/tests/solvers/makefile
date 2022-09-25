BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/pspc/include.mk
include $(SRC_DIR)/pspc/tests/solvers/sources.mk

TEST=pspc/tests/solvers/Test

all: $(pspc_tests_solvers_OBJS) $(BLD_DIR)/$(TEST)

includes:
	@echo $(INCLUDES)

libs:
	@echo $(LIBS)

run: $(pspc_tests_solvers_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pspc/tests/ 

run-quiet: $(pspc_tests_solvers_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pspc/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pspc/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count

clean:
	rm -f $(pspc_tests_solvers_OBJS) $(pspc_tests_solvers_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	$(MAKE) clean-outputs

-include $(pspc_tests_solvers_OBJS:.o=.d)
-include $(pspc_tests_solvers_OBJS:.o=.d)
