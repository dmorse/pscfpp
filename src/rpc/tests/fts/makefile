BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/rpc/include.mk
include $(SRC_DIR)/rpc/tests/fts/sources.mk

TEST=rpc/tests/fts/Test

all: $(rpc_tests_fts_OBJS) $(BLD_DIR)/$(TEST)

includes:
	@echo $(INCLUDES)

libs:
	@echo $(LIBS)

run: $(rpc_tests_fts_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpc/tests/fts/

run-quiet: $(rpc_tests_fts_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpc/tests/fts/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in rpc/tests/fts/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	rm -f out/*

clean:
	rm -f $(BLD_DIR)/$(TEST)
	rm -f $(rpc_tests_fts_OBJS) $(rpc_tests_fts_OBJS:.o=.d)
	$(MAKE) clean-outputs

-include $(rpc_tests_fts_OBJS:.o=.d)
-include $(rpc_tests_fts_OBJS:.o=.d)
