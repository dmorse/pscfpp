
BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/r1d/include.mk
include $(SRC_DIR)/r1d/tests/sources.mk

TEST=r1d/tests/Test

all: $(r1d_tests_OBJS) $(BLD_DIR)/$(TEST)

run: $(r1d_tests_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/r1d/tests/ 

run-quiet: $(r1d_tests_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/r1d/tests/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in r1d/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 
ifeq ($(BLD_DIR),$(SRC_DIR))
	rm -f out/*
else
	rm -f $(SRC_DIR)/r1d/tests/out/*
endif

clean:
	rm -f $(r1d_tests_OBJS) $(r1d_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	$(MAKE) clean-outputs

-include $(r1d_tests_OBJS:.o=.d)
-include $(r1d_OBJS:.o=.d)
