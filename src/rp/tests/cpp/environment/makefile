include ../../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/rp/sources.mk
include $(SRC_DIR)/rp/patterns.mk

.PHONY: all run run-quiet clean-outputs clean

all: Test.o Test

run: Test.o Test
	./Test $(SRC_DIR)/rp/tests/cpp/

run-quiet: Test.o Test
	./Test $(SRC_DIR)/rp/tests/cpp/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rp/tests/cpp/log" > count
	@cat count

clean-outputs:
	rm -f log count out/*.*
	rm -f out/sweep/*

clean:
	rm -f Test
	rm -f Test.o $(rpc_tests_environment_OBJS:.o=.d)
	$(MAKE) clean-outputs

-include $(TEST).d
-include $(rp_DEPS)
-include $(prdc_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
