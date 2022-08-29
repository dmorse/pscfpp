BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/pspg/include.mk
include $(SRC_DIR)/pspg/tests/sources.mk

TEST=pspg/tests/Test

all: $(pspg_tests_OBJS) $(BLD_DIR)/$(TEST)

includes:
	echo $(INCLUDES)

run: $(pspg_tests_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/pspg/tests/ 

run-quiet: $(pspg_tests_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/pspg/tests/  > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pspg/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 
ifeq ($(SRC_DIR),$(BLD_DIR))
	-cd system; $(MAKE) clean-outputs
	-cd field; $(MAKE) clean-outputs
endif

clean:
	rm -f $(pspg_tests_OBJS) $(pspg_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	rm -f log count 
ifeq ($(SRC_DIR),$(BLD_DIR))
	-cd system; $(MAKE) clean
	-cd field; $(MAKE) clean
	rm -f */Test */*.o */*.d 
endif

-include $(pspg_tests_OBJS:.o=.d)
-include $(pspg_OBJS:.o=.d)
