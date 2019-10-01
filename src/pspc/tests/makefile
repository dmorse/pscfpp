BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/pspc/include.mk
include $(SRC_DIR)/pspc/tests/sources.mk

TEST=pspc/tests/Test

all: $(pspc_tests_OBJS) $(BLD_DIR)/$(TEST)

includes:
	echo $(INCLUDES)

run: $(pspc_tests_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pspc/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pspc/tests/log" > count
	@cat count

clean:
	rm -f $(pspc_tests_OBJS) $(pspc_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	rm -f field/Test field/Test.o field/Test.d
	rm -f log count omega.o omegaGrid.o

-include $(pspc_tests_OBJS:.o=.d)
-include $(pspc_OBJS:.o=.d)
