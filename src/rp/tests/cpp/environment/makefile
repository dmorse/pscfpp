include ../../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/rp/sources.mk
include $(SRC_DIR)/rp/patterns.mk

.PHONY: all run run-quiet clean-outputs clean

TEST= $(BLD_DIR)/rp/tests/cpp/environment/Test

all: $(TEST)

$(TEST): $(TEST).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)

run: $(TEST)
	$(TEST) $(SRC_DIR)/rp/tests/cpp/environment/

run-quiet: $(TEST)
	$(TEST) $(SRC_DIR)/rp/tests/cpp/environment > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rp/tests/cpp/environment/log" > count
	@cat count

clean-outputs:
	rm -f log count out/*.*
	rm -f out/sweep/*

clean:
	rm -f $(TEST)
	rm -f $(TEST).o $(TEST).d
	$(MAKE) clean-outputs

-include $(TEST).d
-include $(rp_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
