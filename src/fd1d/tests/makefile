BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/pscf/config.mk
include $(BLD_DIR)/fd1d/config.mk
include $(SRC_DIR)/fd1d/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/fd1d/sources.mk
include $(SRC_DIR)/fd1d/tests/sources.mk

TEST=fd1d/tests/Test

all: $(fd1d_tests_OBJS) $(BLD_DIR)/$(TEST)

run: $(fd1d_tests_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/fd1d/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in fd1d/tests/log" > count
	@cat count

clean:
	rm -f $(fd1d_tests_OBJS) $(fd1d_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	rm -f log count 
ifeq ($(BLD_DIR),$(SRC_DIR))
	rm out/*
endif

-include $(fd1d_tests_OBJS:.o=.d)
-include $(fd1d_OBJS:.o=.d)
