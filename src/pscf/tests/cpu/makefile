include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

TEST=$(BLD_DIR)/pscf/tests/cpu/Test

all: $(TEST).o $(TEST)

$(TEST): $(TEST).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)

run: $(TEST).o $(TEST)
	$(TEST) $(SRC_DIR)/pscf/tests/cpu 

run-quiet: $(TEST).o $(TEST)
	@$(TEST) $(SRC_DIR)/pscf/tests/cpu > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pscf/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 

clean:
	rm -f $(TEST).o
	rm -f $(TEST).d
	rm -f $(TEST)
	$(MAKE) clean-outputs

-include $(TEST).d
-include $(util_DEPS)
-include $(pscf_DEPS)
