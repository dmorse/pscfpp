include ../../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

TEST=$(BLD_DIR)/pscf/tests/backend/cpp/Test

all: $(TEST)

$(TEST): $(TEST).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)

run: $(TEST).o $(TEST)
	$(TEST) $(SRC_DIR)/pscf/tests/backend/cpp/

run-quiet: $(TEST).o $(TEST)
	@$(TEST) $(SRC_DIR)/pscf/tests/backend/cpp/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pscf/tests/cpp/backend/log" > count
	@cat count

clean-outputs:
	rm -f log count 

clean:
	rm -f $(TEST).o
	rm -f $(TEST).d
	rm -f $(TEST)
	$(MAKE) clean-outputs
	rm -f out/*

-include $(TEST).d
-include $(util_DEPS)
-include $(pscf_DEPS)
