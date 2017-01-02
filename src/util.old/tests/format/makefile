BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/util/tests/format/sources.mk

TEST=util/tests/format/TEST

all: $(util_tests_format_OBJS) $(BLD_DIR)/$(TEST)

clean:
	rm -f $(util_tests_format_OBJS) $(util_tests_format_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) 

-include $(util_tests_format_OBJS:.o=.d)

