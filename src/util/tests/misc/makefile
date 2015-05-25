BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/util/tests/misc/sources.mk

TEST=util/tests/misc/Test

all: $(util_tests_misc_OBJS) $(BLD_DIR)/$(TEST)

clean:
	rm -f $(util_tests_misc_OBJS) $(util_tests_misc_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) 

-include $(util_tests_misc_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

