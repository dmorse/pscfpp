include ../../../config.mk
include $(SRC_DIR)/util/config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/tests/cuda/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

TEST=pscf/tests/cuda/Test

all: $(pscf_tests_cuda_OBJS) $(BLD_DIR)/$(TEST)

includes:
	echo $(INCLUDES)

run: $(pscf_tests_cuda_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pscf/tests/cuda 

run-quiet: $(pscf_tests_cuda_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/pscf/tests/cuda > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pscf/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 

clean:
	rm -f $(pscf_tests_cuda_OBJS) $(pscf_tests_cuda_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST)
	$(MAKE) clean-outputs

-include $(pscf_tests_cuda_OBJS:.o=.d)
-include $(pscf_util_OBJS:.o=.d)
-include $(pscf_pscf_OBJS:.o=.d)
