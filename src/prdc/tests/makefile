include ../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/prdc/patterns.mk

TEST_CUDA=$(BLD_DIR)/prdc/tests/cudaTest
TEST_CPU=$(BLD_DIR)/prdc/tests/cpuTest

ifdef PSCF_CUDA
TEST=$(TEST_CUDA)
TEST_OBJ=$(TEST).ou
TEST_DEP=$(TEST).du
else
TEST=$(TEST_CPU)
TEST_OBJ=$(TEST).o
TEST_DEP=$(TEST).d
endif

all: $(TEST) $(PSCF_LIB)

$(TEST_CPU): $(TEST_CPU).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(TEST_CPU) $(TEST_CPU).o $(LIBS)

$(TEST_CUDA): $(TEST_CUDA).ou $(PSCF_LIBS)
	@cp $(TEST_CUDA).ou $(TEST_CUDA).o
	$(NVXX) $(LDFLAGS) -o $(TEST_CUDA) $(TEST_CUDA).o $(LIBS)
	@rm $(TEST_CUDA).o

lib:
	echo $(PSCF_LIBS)

run: $(TEST)
	@$(TEST) $(SRC_DIR)/prdc/tests/ 

run-quiet: $(TEST_OBJ) $(TEST)
	@$(TEST) $(SRC_DIR)/prdc/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in prdc/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	-cd $(SRC_DIR)/prdc/tests/crystal; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/prdc/tests/cpu; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/prdc/tests/cuda; $(MAKE) clean-outputs

clean:
	rm -f $(TEST_OBJ)
	rm -f $(TEST_DEP)
	rm -f $(TEST)
	rm -f cpuTest cudaTest
	rm -f *Test.o *Test.d
	rm -f *Test.ou *Test.du
	rm -f log count
	-cd $(SRC_DIR)/prdc/tests/crystal; $(MAKE) clean
	-cd $(SRC_DIR)/prdc/tests/cpu; $(MAKE) clean
	-cd $(SRC_DIR)/prdc/tests/cuda; $(MAKE) clean

.PHONY: all run run-quiet clean-outputs clean

-include $(util_DEPS)
-include $(pscf_DEPS)
-include $(prdc_DEPS)
-include $(TEST_DEP)
