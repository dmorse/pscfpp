BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/cyln/include.mk
include $(SRC_DIR)/cyln/tests/field/sources.mk

TEST=cyln/tests/field/Test

all: $(cyln_tests_field_OBJS) $(BLD_DIR)/$(TEST)

includes:
	@echo $(INCLUDES)

libs:
	@echo $(LIBS)

run: $(cyln_tests_field_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/cyln/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in cyln/tests/log" > count
	@cat count

clean:
	rm -f $(cyln_tests_field_OBJS) $(cyln_tests_field_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	rm -f log count binary

-include $(cyln_tests_field_OBJS:.o=.d)
-include $(cyln_tests_field_OBJS:.o=.d)
