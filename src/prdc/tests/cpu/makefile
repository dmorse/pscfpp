BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/prdc/include.mk
include $(SRC_DIR)/prdc/tests/cpu/sources.mk

TEST=prdc/tests/cpu/Test

all: $(prdc_tests_cpu_OBJS) $(BLD_DIR)/$(TEST)

includes:
	@echo $(INCLUDES)

libs:
	@echo $(LIBS)

run: $(prdc_tests_cpu_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/prdc/tests/ 

run-quiet: $(prdc_tests_cpu_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/prdc/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in prdc/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count out/*

clean:
	rm -f $(prdc_tests_cpu_OBJS) $(prdc_tests_cpu_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	$(MAKE) clean-outputs

-include $(prdc_tests_cpu_OBJS:.o=.d)
-include $(prdc_tests_cpu_OBJS:.o=.d)
