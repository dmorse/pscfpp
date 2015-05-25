BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/util/tests/boundary/sources.mk

TEST=util/tests/boundary/Test

all: $(util_tests_boundary_OBJS) $(BLD_DIR)/$(TEST)

clean:
	rm -f $(util_tests_boundary_OBJS) $(BLD_DIR)/$(TEST) $(util_tests_boundary_OBJS:.o=.d)

clean-deps:
	rm -f $(util_tests_boundary_OBJS:.o=.d)

-include $(util_tests_boundary_OBJS:.o=.d)

