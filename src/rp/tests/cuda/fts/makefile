include ../../../../config.mk
include $(SRC_DIR)/rp/include.mk

TEST=$(BLD_DIR)/rp/tests/cuda/fts/Test

all: $(TEST).ou $(TEST)

$(TEST): $(TEST).ou $(PSCF_LIBS)
	cp $(TEST).ou $(TEST).o
	$(NVXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)
	rm $(TEST).o

run: $(TEST).ou $(TEST)
	$(TEST) $(SRC_DIR)/rpg/tests/fts/

run-quiet: $(TEST).ou $(TEST)
	$(TEST) $(SRC_DIR)/rpg/tests/fts/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in rpg/tests/fts/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	rm -f out/*

clean:
	rm -f $(TEST)
	rm -f $(TEST).ou 
	rm -f $(TEST).ou 
	$(MAKE) clean-outputs

-include $(TEST).du
