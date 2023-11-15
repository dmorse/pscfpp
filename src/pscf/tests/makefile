BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/pscf/include.mk
include $(SRC_DIR)/pscf/tests/sources.mk

ifdef PSCF_CUDA
TEST=pscf/tests/cudaTest
else
TEST=pscf/tests/cpuTest
endif

.PHONY: all libs includes run run-quiet clean-outputs clean

all: $(pscf_tests_OBJS) $(BLD_DIR)/$(TEST)

libs:
	@echo $(LIBS)

includes:
	@echo $(INCLUDES)

run: $(pscf_tests_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pscf/tests/ 

run-quiet: $(pscf_tests_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/pscf/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pscf/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 
ifeq ($(BLD_DIR),$(SRC_DIR))
	-cd chem; $(MAKE) clean-outputs
	-cd inter; $(MAKE) clean-outputs
	-cd math; $(MAKE) clean-outputs
	-cd mesh; $(MAKE) clean-outputs
	-cd solvers; $(MAKE) clean-outputs
else
	-cd $(SRC_DIR)/pscf/tests/chem; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/pscf/tests/inter; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/pscf/tests/math; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/pscf/tests/mesh; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/pscf/tests/solvers; $(MAKE) clean-outputs
endif

clean:
	rm -f *Test
	rm -f $(pscf_tests_OBJS) $(pscf_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
ifeq ($(BLD_DIR),$(SRC_DIR))
	rm -f chem/Test chem/Test.o chem/Test.d
	rm -f inter/Test inter/Test.o inter/Test.d
	rm -f math/Test math/Test.o math/Test.d
	rm -f mesh/Test mesh/Test.o mesh/Test.d
	rm -f solvers/Test solvers/Test.o solvers/Test.d
else
	-cd $(SRC_DIR)/pscf/tests/chem; $(MAKE) clean
	-cd $(SRC_DIR)/pscf/tests/inter; $(MAKE) clean
	-cd $(SRC_DIR)/pscf/tests/math; $(MAKE) clean
	-cd $(SRC_DIR)/pscf/tests/mesh; $(MAKE) clean
	-cd $(SRC_DIR)/pscf/tests/solvers; $(MAKE) clean
endif
	$(MAKE) clean-outputs

-include $(pscf_tests_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
