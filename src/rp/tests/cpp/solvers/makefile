# -------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/rp/sources.mk
include $(SRC_DIR)/rp/patterns.mk

# -------------------------------------------------------------------------
# Main targets

TEST=$(BLD_DIR)/rp/tests/cpp/solvers/Test

all: $(TEST)

$(TEST): $(TEST).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)

run: $(TEST).o $(TEST)
	./$(TEST) $(SRC_DIR)/rp/tests/cpp/solvers 

run-quiet: $(TEST)
	@./$(TEST) $(SRC_DIR)/rp/tests/cpp/solvers/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rp/tests/cpp/solvers/log" > count
	@cat count

clean-outputs:
	rm -f log count

clean:
	rm -f $(TEST)
	rm -f $(TEST).o $(TEST).d
	$(MAKE) clean-outputs

.PHONY= all run run-quiet clean-outputs clean
# -------------------------------------------------------------------------
# Main targets

-include Test.d
-include $(rp_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
