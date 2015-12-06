BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/pfts/config.mk
include $(SRC_DIR)/pfts/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pfts/sources.mk
include $(SRC_DIR)/pfts/tests/solvers/sources.mk

TEST=pfts/tests/solvers/Test

all: $(pfts_tests_solvers_OBJS) $(BLD_DIR)/$(TEST)

clean:
	rm -f $(pfts_tests_solvers_OBJS) $(pfts_tests_solvers_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) 

-include $(pfts_tests_solvers_OBJS:.o=.d)
-include $(pfts_OBJS:.o=.d)

