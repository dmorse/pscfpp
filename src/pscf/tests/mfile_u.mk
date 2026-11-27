
# This file is included by makefiles in subdirectories of the directory
# src/pscf/tests that contain a CUDA C++ file named Test.cu
#
# Requirements:
#  - The src/config.mk file must be included before this file
#  - Variable SHORT_DIR must be defined before this file is included
#  - Target clean-outputs must be defined after this file is included

# -------------------------------------------------------------------
# Include makefile fragments
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk


# -------------------------------------------------------------------
# Definitions

TEST_DIR=$(SRC_DIR)/$(SHORT_DIR)
TEST=$(BLD_DIR)/$(SHORT_DIR)/Test

# -------------------------------------------------------------------
# Main targets

.PHONY= all run run-quiet clean clean-outputs

all: $(TEST)

$(TEST): $(TEST).ou $(PSCF_LIBS)
	@cp $(TEST).ou $(TEST).o
	$(NVXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)
	@rm $(TEST).o

run: $(TEST)
	$(TEST) $(TEST_DIR)/

run-quiet: $(TEST)
	@$(TEST) $(TEST_DIR)/  > log
	@echo `grep failed log` ", "\
              `grep successful log` "in $(SHORT_DIR)/log" > count
	@cat count

clean:
	rm -f $(TEST)
	rm -f $(TEST).ou
	rm -f $(TEST).du
	$(MAKE) clean-outputs

# -------------------------------------------------------------------
# Including dependency files

-include $(TEST).du
-include $(pscf_DEPS)
-include $(util_DEPS)
