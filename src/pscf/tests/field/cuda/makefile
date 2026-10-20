include ../../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/tests/field/cuda/sources.mk

TEST=$(BLD_DIR)/pscf/tests/field/cuda/Test

all: $(TEST)

$(TEST): $(TEST).ou $(PSCF_LIBS)
	@cp $(TEST).ou $(TEST).o
	$(NVXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)
	@rm $(TEST).o

run: $(TEST)
	$(TEST) $(SRC_DIR)/pscf/tests/field/cuda/

run-quiet: $(TEST)
	@$(TEST) $(SRC_DIR)/pscf/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pscf/tests/field/cuda/log" > count
	@cat count

clean-outputs:
	rm -f log count 

clean:
	rm -f $(TEST) $(TEST).ou $(TEST).du
	$(MAKE) clean-outputs

.PHONY: all run run-quiet clean-outputs clean

-include $(pscf_DEPS)
-include $(util_DEPS)
-include $(TEST).du
