include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/rp/sources.mk
include $(SRC_DIR)/rp/patterns.mk

# Definitions

# -------------------------------------------------------------------
# Definitions

SHORT_DIR=rp/tests/cuda
TEST_DIR=$(SRC_DIR)/$(SHORT_DIR)
TEST=$(BLD_DIR)/$(SHORT_DIR)/Test

# -------------------------------------------------------------------
# Main targets

.PHONY= all run run-quiet clean clean-outputs

all: $(TEST)

$(TEST): $(TEST).ou $(PSCF_LIBS)
	@cp $(TEST).ou $(TEST).o
	$(NVXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)
	@rm $(TEST).o

run: $(TEST)
	$(TEST) $(TEST_DIR)/

run-quiet: $(TEST)
	@$(TEST) $(TEST_DIR)/  > log
	@echo `grep failed log` ", "\
              `grep successful log` "in $(SHORT_DIR)/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	-cd $(TEST_DIR)/field; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/environment; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/solvers; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/system; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/sweep; $(MAKE) clean-outputs
	-cd $(TEST_DIR)/fts; $(MAKE) clean-outputs

clean:
	rm -f $(TEST)
	rm -f $(TEST).ou
	rm -f $(TEST).du
	rm -f log count 
	-cd $(TEST_DIR)/field; $(MAKE) clean
	-cd $(TEST_DIR)/environment; $(MAKE) clean
	-cd $(TEST_DIR)/solvers; $(MAKE) clean
	-cd $(TEST_DIR)/system; $(MAKE) clean
	-cd $(TEST_DIR)/sweep; $(MAKE) clean
	-cd $(TEST_DIR)/fts; $(MAKE) clean
	rm -f $(TEST_DIR)/*/Test
	rm -f $(TEST_DIR)/*.ou $(TEST_DIR)/*.du
	rm -f $(TEST_DIR)/*/*.ou $(TEST_DIR)/*/*.du

.PHONY: all run run-quiet clean-outputs clean

-include $(TEST).du
-include $(rp_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
