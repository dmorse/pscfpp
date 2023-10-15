BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/prdc/include.mk
include $(SRC_DIR)/prdc/tests/sources.mk

.PHONY: all objs list includes run run-quiet clean-outputs clean

TEST=prdc/tests/Test

all: $(prdc_tests_OBJS) $(BLD_DIR)/$(TEST)

list:
	@echo $(prdc_tests_)
	@echo $(prdc_tests_OBJS)

includes:
	@echo $(INCLUDES)

run: $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/prdc/tests/ 

run-quiet: $(prdc_tests_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/prdc/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in prdc/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 
ifeq ($(BLD_DIR),$(SRC_DIR))
	-cd crystal; $(MAKE) clean-outputs
else
	-cd $(SRC_DIR)/prdc/tests/crystal; $(MAKE) clean-outputs
endif

clean:
	rm -f $(prdc_tests_OBJS) $(prdc_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	rm -f Test.cpp
ifeq ($(BLD_DIR),$(SRC_DIR))
	rm -f crystal/Test crystal/Test.o crystal/Test.d
else
	-cd $(SRC_DIR)/prdc/tests/crystal; $(MAKE) clean
endif
	$(MAKE) clean-outputs

-include $(prdc_tests_OBJS:.o=.d)
-include $(prdc_OBJS:.o=.d)
