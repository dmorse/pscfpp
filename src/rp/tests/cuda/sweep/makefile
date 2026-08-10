include ../../../../config.mk
include $(SRC_DIR)/rp/include.mk

TEST=$(BLD_DIR)/rp/tests/cuda/sweep/Test

all: $(TEST)

$(TEST): $(TEST).ou $(PSCF_LIBS)
	cp $(TEST).ou $(TEST).o
	$(NVXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)
	rm $(TEST).o

run: $(TEST).ou $(TEST)
	@$(TEST) $(SRC_DIR)/rp/tests/cuda

run-quiet: $(TEST).ou $(TEST)
	@$(TEST) $(SRC_DIR)/rp/tests/cuda  > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rpg/tests/sweep/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	find out -type f \( ! -iname ".*" \) -exec rm -f {} \;

clean:
	rm -f $(TEST)
	rm -f $(TEST).ou 
	rm -f $(TEST).du 
	$(MAKE) clean-outputs

-include $(TEST).du
