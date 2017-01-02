BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/util/tests/random/sources.mk

TEST=util/tests/random/Test

all: $(util_tests_random_OBJS) $(BLD_DIR)/$(TEST)

clean:
	rm -f $(util_tests_random_OBJS) $(util_tests_random_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) 
	$(MAKE) clean-outputs

clean-outputs:
	rm -f binary

-include $(util_tests_random_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

