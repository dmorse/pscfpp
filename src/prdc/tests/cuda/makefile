include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/tests/cuda/sources.mk
include $(SRC_DIR)/prdc/patterns.mk

TEST=$(BLD_DIR)/prdc/tests/cuda/Test

all: $(prdc_tests_cuda_OBJS) $(TEST)

$(TEST): $(TEST).ou $(PSCF_LIBS)
	@cp $(TEST).ou $(TEST).o
	$(NVXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)
	@rm $(TEST).o

run: $(prdc_tests_cuda_OBJS) $(TEST)
	$(TEST) $(SRC_DIR)/prdc/tests/cuda/

run-quiet: $(prdc_tests_cuda_OBJS) $(TEST)
	@$(TEST) $(SRC_DIR)/prdc/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in prdc/tests/cuda/log" > count
	@cat count

clean-outputs:
	rm -f log count 

clean:
	rm -f $(prdc_tests_cuda_OBJS)
	rm -f $(prdc_tests_cuda_DEPS)
	rm -f $(TEST)
	$(MAKE) clean-outputs

.PHONY: all run run-quiet clean-outputs clean

-include $(prdc_tests_cuda_DEPS)
-include $(prdc_DEPS)
-include $(util_DEPS)
-include $(pscf_DEPS)
-include $(TEST).du
