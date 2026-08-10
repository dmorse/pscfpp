include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

TEST=$(BLD_DIR)/pscf/tests/interaction/Test

all: $(TEST)

run: $(TEST).o $(TEST)
	$(TEST) $(SRC_DIR)/pscf/tests/interaction/ 

run-quiet: $(TEST)
	@$(TEST) $(SRC_DIR)/pscf/tests/interaction/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pscf/tests/interaction/log" > count
	@cat count

clean-outputs:
	rm -f log count 

clean:
	rm -f $(TEST).o 
	rm -f $(TEST).d
	rm -f $(TEST)
	$(MAKE) clean-outputs

-include $(TEST).d
-include $(pscf_DEPS)
-include $(util_DEPS)
