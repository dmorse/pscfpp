include ../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

TEST_CUDA=$(BLD_DIR)/pscf/tests/cudaTest
TEST_CPU=$(BLD_DIR)/pscf/tests/cpuTest
TEST_DIR=$(SRC_DIR)/pscf/tests

ifdef PSCF_CUDA
TEST=$(TEST_CUDA)
TEST_OBJ=$(TEST).ou
TEST_DEP=$(TEST).du
else
TEST=$(TEST_CPU)
TEST_OBJ=$(TEST).o
TEST_DEP=$(TEST).d
endif

all: $(TEST)

$(TEST_CPU): $(TEST_CPU).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(TEST_CPU) $(TEST_CPU).o $(LIBS)

$(TEST_CUDA): $(TEST_CUDA).ou $(PSCF_LIBS)
	@cp $(TEST_CUDA).ou $(TEST_CUDA).o
	$(NVXX) $(LDFLAGS) -o $(TEST_CUDA) $(TEST_CUDA).o $(LIBS)
	@rm $(TEST_CUDA).o

run: $(TEST)
	@$(TEST) $(SRC_DIR)/pscf/tests/ 

run-quiet: $(TEST_OBJ) $(TEST)
	@$(TEST) $(SRC_DIR)/pscf/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pscf/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	-cd $(TEST_DIR)/chem; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/floryHuggins; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/correlation; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/interaction; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/math; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/mesh; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/solvers; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/cpu; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/cuda; $(MAKE) clean-outputs

clean:
	rm -f $(TEST_OBJ)
	rm -f $(TEST_DEP)
	rm -f $(TEST)
	rm -f cpuTest cudaTest
	rm -f *.o *.d 
	rm -f *.ou *.du
	rm -f log count 
	-cd $(TEST_DIR)/chem; $(MAKE) clean
	-cd $(TEST_DIR)/floryHuggins; $(MAKE) clean
	-cd $(TEST_DIR)/correlation; $(MAKE) clean
	-cd $(TEST_DIR)/interaction; $(MAKE) clean
	-cd $(TEST_DIR)/math; $(MAKE) clean
	-cd $(TEST_DIR)/mesh; $(MAKE) clean
	-cd $(TEST_DIR)/solvers; $(MAKE) clean
	-cd $(TEST_DIR)/cpu; $(MAKE) clean
	-cd $(TEST_DIR)/cuda; $(MAKE) clean

.PHONY: all run run-quiet clean-outputs clean

# Include all relevant dependency files
-include $(util_DEPS)
-include $(pscf_DEPS)
-include $(TEST_DEP)
