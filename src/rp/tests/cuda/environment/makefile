include ../../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/rp/sources.mk
include $(SRC_DIR)/rp/patterns.mk

TEST=$(BLD_DIR)/rp/tests/cuda/environment/Test

all: $(TEST)

$(TEST): $(TEST).ou $(PSCF_LIBS)
	@cp $(TEST).ou $(TEST).o
	$(NVXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)
	@rm $(TEST).o

run: $(TEST).ou $(TEST)
	$(TEST) $(SRC_DIR)/rp/tests/cuda/ 

run-quiet: $(TEST).ou $(TEST)
	$(TEST) $(SRC_DIR)/rp/tests/cuda/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rp/tests/cuda/log" > count
	@cat count

clean-outputs:
	rm -f log count out/*.*
	rm -f out/sweep/*

clean:
	rm -f $(TEST)
	rm -f $(TEST).ou 
	rm -f $(TEST).du 
	$(MAKE) clean-outputs

.PHONY: all run run-quiet clean-outputs clean

-include $(TEST).du
-include $(rp_DEPS)
-include $(prdc_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
