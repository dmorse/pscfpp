BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/prdc/include.mk
include $(SRC_DIR)/prdc/tests/sources.mk

.PHONY: all objs list includes run run-quiet clean-outputs clean

ifdef PSCF_CUDA
TEST=prdc/tests/cudaTest
else
TEST=prdc/tests/cpuTest
endif

all: $(prdc_tests_OBJS) $(BLD_DIR)/$(TEST)

list:
	@echo $(prdc_tests_)
	@echo $(prdc_tests_OBJS)

includes:
	@echo $(INCLUDES)

run: $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/prdc/tests/ 

run-quiet: $(prdc_tests_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/prdc/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in prdc/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	-cd $(SRC_DIR)/prdc/tests/crystal; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/prdc/tests/cpu; $(MAKE) clean-outputs
	-cd $(SRC_DIR)/prdc/tests/cuda; $(MAKE) clean-outputs

clean:
	rm -f $(prdc_tests_OBJS) $(prdc_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST)
	rm -f cpuTest cudaTest
	rm -f *Test.o *Test.d
	rm -f log count
	-cd $(SRC_DIR)/prdc/tests/crystal; $(MAKE) clean
	-cd $(SRC_DIR)/prdc/tests/cpu; $(MAKE) clean
	-cd $(SRC_DIR)/prdc/tests/cuda; $(MAKE) clean

-include $(prdc_tests_OBJS:.o=.d)
-include $(prdc_OBJS:.o=.d)
