BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/chem/config.mk
include $(SRC_DIR)/chem/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/chem/sources.mk
include $(SRC_DIR)/chem/tests/sources.mk

TEST=chem/tests/Test

all: $(chem_tests_OBJS) $(BLD_DIR)/$(TEST)

run: $(chem_tests_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/chem/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in chem/tests/log" > count
	@cat count

clean:
	rm -f $(chem_tests_OBJS) $(chem_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	rm -f log count 

-include $(chem_tests_OBJS:.o=.d)
-include $(chem_OBJS:.o=.d)
