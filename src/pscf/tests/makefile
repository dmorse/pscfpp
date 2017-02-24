BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/pscf/config.mk
include $(SRC_DIR)/pscf/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
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
	rm -f math/Test math/Test.o math/Test.d
	rm -f chem/Test chem/Test.o chem/Test.d
	rm -f solvers/Test solvers/Test.o solvers/Test.d
	rm -f inter/Test inter/Test.o inter/Test.d
	rm -f mesh/Test mesh/Test.o mesh/Test.d
	rm -f crystal/Test crystal/Test.o crystal/Test.d
	rm -f log count 

-include $(pscf_tests_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
