include ../../../config.mk
include $(SRC_DIR)/util/config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/tests/cuda/sources.mk
include $(SRC_DIR)/prdc/patterns.mk

TEST=prdc/tests/cuda/Test

all: $(prdc_tests_cuda_OBJS) $(BLD_DIR)/$(TEST)

run: $(prdc_tests_cuda_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/prdc/tests/cuda/

run-quiet: $(prdc_tests_cuda_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/prdc/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in prdc/tests/cuda/log" > count
	@cat count

clean-outputs:
	rm -f log count 

clean:
	rm -f $(prdc_tests_cuda_OBJS) $(prdc_tests_cuda_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	$(MAKE) clean-outputs

-include $(prdc_tests_cuda_OBJS:.o=.d)
-include $(prdc_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
