# ---------------------------------------------------------------------
# Include makefile fragments

include ../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

# ---------------------------------------------------------------------
# Variable definitions

SHORT_DIR=pscf/tests

TEST_DIR=$(SRC_DIR)/$(SHORT_DIR)
TEST=$(BLD_DIR)/pscf/tests/Test

ifdef PSCF_CUDA
TEST_OBJ=$(TEST).ou
TEST_DEP=$(TEST).du
else
TEST_OBJ=$(TEST).o
TEST_DEP=$(TEST).d
endif

# ---------------------------------------------------------------------
# Main targets

.PHONY: all run run-quiet clean clean-outputs

all: $(TEST)

ifdef PSCF_CUDA

$(TEST): $(TEST).ou $(PSCF_LIBS)
	@cp $(TEST).ou $(TEST).o
	$(NVXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)
	@rm $(TEST).o

else

$(TEST): $(TEST).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)

endif

run: $(TEST)
	@$(TEST) $(TEST_DIR)/ 

run-quiet: $(TEST_OBJ) $(TEST)
	@$(TEST) $(TEST_DIR)/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in $(SHORT_DIR)" > count
	@cat count

clean:
	rm -f $(TEST_OBJ)
	rm -f $(TEST_DEP)
	rm -f $(TEST)
	rm -f cppTest cudaTest
	rm -f *.o *.d
	rm -f *.ou *.du
	rm -f log count 
	-cd $(TEST_DIR)/backend/cpp; $(MAKE) clean
	-cd $(TEST_DIR)/backend/cuda; $(MAKE) clean
	-cd $(TEST_DIR)/chem; $(MAKE) clean
	-cd $(TEST_DIR)/correlation; $(MAKE) clean
	-cd $(TEST_DIR)/crystal; $(MAKE) clean
	-cd $(TEST_DIR)/floryHuggins; $(MAKE) clean
	-cd $(TEST_DIR)/interaction; $(MAKE) clean
	-cd $(TEST_DIR)/math; $(MAKE) clean
	-cd $(TEST_DIR)/mesh; $(MAKE) clean
	-cd $(TEST_DIR)/solvers; $(MAKE) clean

clean-outputs:
	rm -f log count 
	-cd $(TEST_DIR)/backend/cpp; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/backend/cuda; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/chem; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/correlation; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/crystal; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/floryHuggins; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/interaction; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/math; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/mesh; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/solvers; $(MAKE) clean-outputs

# ---------------------------------------------------------------------
# Include dependency files

-include $(TEST_DEP)
-include $(util_DEPS)
-include $(pscf_DEPS)
