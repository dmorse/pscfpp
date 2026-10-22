
include ../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/r1/sources.mk
include $(SRC_DIR)/r1/patterns.mk

TEST=$(BLD_DIR)/r1/tests/Test

all: $(TEST)

$(TEST): $(TEST).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)

run: $(TEST)
	$(TEST) $(SRC_DIR)/r1/tests/ 

run-quiet: $(TEST)
	@$(TEST) $(SRC_DIR)/r1/tests/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in r1/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	rm -f $(SRC_DIR)/r1/tests/out/*

clean:
	rm -f $(TEST) 
	rm -f $(TEST).o $(TEST).d
	$(MAKE) clean-outputs

-include $(util_DEPS)
-include $(pscf_DEPS)
-include $(r1_DEPS)
-include $(TEST).d
