BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/pscf/include.mk
include $(SRC_DIR)/pscf/tests/sources.mk

TEST=pscf/tests/Test

all: $(pscf_tests_OBJS) $(BLD_DIR)/$(TEST)

includes:
	echo $(INCLUDES)

run: $(pscf_tests_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pscf/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pscf/tests/log" > count
	@cat count

clean:
	rm -f $(pscf_tests_OBJS) $(pscf_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	rm -f chem/Test chem/Test.o chem/Test.d
	rm -f crystal/Test crystal/Test.o crystal/Test.d
	rm -f inter/Test inter/Test.o inter/Test.d
	rm -f math/Test math/Test.o math/Test.d
	rm -f mesh/Test mesh/Test.o mesh/Test.d
	rm -f solvers/Test solvers/Test.o solvers/Test.d
	rm -f log count 

clean-outputs:
	rm -f log count 
	cd chem; $(MAKE) clean-outputs
	cd crystal; $(MAKE) clean-outputs
	cd inter; $(MAKE) clean-outputs
	cd math; $(MAKE) clean-outputs
	cd mesh; $(MAKE) clean-outputs
	cd solvers; $(MAKE) clean-outputs

-include $(pscf_tests_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
