BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/pfts/config.mk
include $(SRC_DIR)/pfts/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pfts/sources.mk
include $(SRC_DIR)/pfts/tests/sources.mk

TEST=pfts/tests/Test

all: $(pfts_tests_OBJS) $(BLD_DIR)/$(TEST)

run: $(pfts_tests_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pfts/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pfts/tests/log" > count
	@cat count

clean:
	rm -f $(pfts_tests_OBJS) $(pfts_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	rm -f log count dummy
ifeq ($(BLD_DIR),$(SRC_DIR))
	cd chem; $(MAKE) clean
	rm -f */*.o */*/*.o 
	rm -f */*.d */*/*.d
else
	cd $(SRC_DIR)/pfts/tests; $(MAKE) clean-outputs
endif

ifeq ($(BLD_DIR),$(SRC_DIR))
clean-outputs:
	cd chem; $(MAKE) clean-outputs
endif

-include $(pfts_tests_OBJS:.o=.d)
-include $(pfts_OBJS:.o=.d)
