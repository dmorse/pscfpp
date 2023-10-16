BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/pspc/include.mk
include $(SRC_DIR)/pspc/tests/field/sources.mk

TEST=pspc/tests/field/Test

.PHONY: all includes defines libs run run-quiet clean-outputs clean

all: $(pspc_tests_field_OBJS) $(BLD_DIR)/$(TEST)

includes:
	@echo $(INCLUDES)

defines:
	@echo $(DEFINES)

libs:
	@echo $(LIBS)

run: $(pspc_tests_field_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pspc/tests/ 

run-quiet: $(pspc_tests_field_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pspc/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pspc/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count out/*

clean:
	rm -f $(pspc_tests_field_OBJS) $(pspc_tests_field_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	$(MAKE) clean-outputs

-include $(pspc_OBJS:.o=.d)
-include $(pspc_tests_field_OBJS:.o=.d)
