# --------------------------------------------------------------------
# Include makefile fragments

include ../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/r1/sources.mk
include $(SRC_DIR)/r1/patterns.mk

# --------------------------------------------------------------------
# Definitions

SHORT_DIR=r1/tests
TEST_DIR=$(SRC_DIR)/$(SHORT_DIR)
TEST=$(BLD_DIR)/$(SHORT_DIR)/Test

# --------------------------------------------------------------------
# Targets and rules

.PHONY = all run run-quiet clean clean-outputs

all: $(TEST)

$(TEST): $(TEST).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)

run: $(TEST)
	$(TEST) $(TEST_DIR)/ 

run-quiet: $(TEST)
	@$(TEST) $(TEST_DIR)/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in $(SHORT_DIR)/log" > count
	@cat count

clean:
	rm -f $(TEST) 
	rm -f $(TEST).o $(TEST).d
	$(MAKE) clean-outputs

clean-outputs:
	rm -f log count 
	rm -f $(TEST_DIR)/out/*

# --------------------------------------------------------------------
# Include dependency files

-include $(util_DEPS)
-include $(pscf_DEPS)
-include $(r1_DEPS)
-include $(TEST).d
