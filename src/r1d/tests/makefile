
include ../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/r1d/sources.mk
include $(SRC_DIR)/r1d/patterns.mk

TEST=$(BLD_DIR)/r1d/tests/Test

all: $(TEST)

$(TEST): $(TEST).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)

run: $(TEST)
	$(TEST) $(SRC_DIR)/r1d/tests/ 

run-quiet: $(TEST)
	@$(TEST) $(SRC_DIR)/r1d/tests/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in r1d/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	rm -f $(SRC_DIR)/r1d/tests/out/*

clean:
	rm -f $(TEST) 
	rm -f $(TEST).o $(TEST).d
	$(MAKE) clean-outputs

-include $(util_DEPS)
-include $(pscf_DEPS)
-include $(r1d_DEPS)
-include $(TEST).d
