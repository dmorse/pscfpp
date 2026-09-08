include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/prdc/patterns.mk

TEST=$(BLD_DIR)/prdc/tests/crystal/Test

all: $(TEST).o $(TEST)

$(TEST): $(TEST).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)

run:    $(TEST).o $(TEST)
	@$(TEST) $(SRC_DIR)/prdc/tests/crystal/ 

run-quiet: $(TEST).o $(TEST)
	@$(TEST) $(SRC_DIR)/prdc/tests/crystal/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in prdc/tests/crystal/log" > count
	@cat count

clean-outputs:
	rm -f log count 

clean:
	rm -f $(TEST).o $(TEST).d
	rm -f $(TEST)
	$(MAKE) clean-outputs

.PHONY= all run run-quiet clean-outputs clean

-include $(TEST).d
-include $(prdc_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
