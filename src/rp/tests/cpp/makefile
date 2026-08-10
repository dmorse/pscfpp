# ------------------------------------------------------------------
#  Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/rp/sources.mk
include $(SRC_DIR)/rp/patterns.mk

TEST=$(BLD_DIR)/rp/tests/cpp/Test

all: $(TEST)

$(TEST): $(TEST).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)

run: $(TEST)
	$(TEST) $(SRC_DIR)/rp/tests/cpp/

run-quiet: $(TEST)
	$(TEST) $(SRC_DIR)/rp/tests/cpp/  > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rp/tests/cpp/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	-cd $(TEST_CPP_DIR)/field; $(MAKE) clean-outputs
	-cd $(TEST_CPP_DIR)/environment; $(MAKE) clean-outputs
	-cd $(TEST_CPP_DIR)/solvers; $(MAKE) clean-outputs
	-cd $(TEST_CPP_DIR)/system; $(MAKE) clean-outputs
	-cd $(TEST_CPP_DIR)/sweep; $(MAKE) clean-outputs
	-cd $(TEST_CPP_DIR)/fts; $(MAKE) clean-outputs

clean:
	rm -f $(TEST) 
	rm -f $(TEST).o $(TEST).d
	rm -f log count 
	-cd $(TEST_CPP_DIR)/field; $(MAKE) clean
	-cd $(TEST_CPP_DIR)/environment; $(MAKE) clean
	-cd $(TEST_CPP_DIR)/solvers; $(MAKE) clean
	-cd $(TEST_CPP_DIR)/system; $(MAKE) clean
	-cd $(TEST_CPP_DIR)/sweep; $(MAKE) clean
	-cd $(TEST_CPP_DIR)/fts; $(MAKE) clean
	rm -f $(TEST_CPP_DIR)/*/Test
	rm -f $(TEST_CPP_DIR)/*.o $(TEST_CPP_DIR)/*.d
	rm -f $(TEST_CPP_DIR)/*/*.o $(TEST_CPP_DIR)/*/*.d

.PHONY= all run run-quiet clean-outputs clean

# ------------------------------------------------------------------
#  Include dependency files

-include $(TEST).d
-include $(rp_DEPS)
-include $(prdc_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
