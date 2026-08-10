include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

TEST=$(BLD_DIR)/pscf/tests/cuda/cudaTest

all: $(TEST)

$(TEST): $(TEST).ou $(PSCF_LIBS)
	@cp $(TEST).ou $(TEST).o
	$(NVXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)
	@rm $(TEST).o

run: $(TEST)
	$(TEST) $(SRC_DIR)/pscf/tests/cuda 

run-quiet: $(TEST)
	@$(TEST) $(SRC_DIR)/pscf/tests/cuda > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pscf/tests/log" > count
	@cat count

clean-outputs:
	rm -f log count 

clean:
	rm -f $(TEST).ou
	rm -f $(TEST).du
	rm -f $(TEST)
	$(MAKE) clean-outputs

-include $(TEST).du
-include $(util_DEPS)
-include $(pscf_DEPS)
