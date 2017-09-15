BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/pscf/config.mk
include $(SRC_DIR)/pscf/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/tests/math/sources.mk

TEST=pscf/tests/math/Test

all: $(pscf_tests_math_OBJS) $(BLD_DIR)/$(TEST)

includes:
	echo $(INCLUDES)

run: $(pscf_tests_math_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pscf/tests/math > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pscf/tests/log" > count
	@cat count

clean:
	rm -f $(pscf_tests_math_OBJS) $(pscf_tests_math_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	rm -f log count 

-include $(pscf_tests_math_OBJS:.o=.d)
-include $(pscf_tests_math_OBJS:.o=.d)
