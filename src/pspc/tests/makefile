BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/pspc/include.mk
include $(SRC_DIR)/pspc/tests/sources.mk

TEST=pspc/tests/Test

all: $(pspc_tests_OBJS) $(BLD_DIR)/$(TEST)

includes:
	echo $(INCLUDES)

run: $(pspc_tests_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pspc/tests/ 

#	@echo `grep failed log` ", "\
#              `grep successful log` "in pspc/tests/log" > count
#	@cat count

clean-outputs:
	rm -f log count 
ifeq ($(SRC_DIR),$(BLD_DIR))
	-cd field; $(MAKE) clean-outputs
	-cd solvers; $(MAKE) clean-outputs
	-cd system; $(MAKE) clean-outputs
endif

clean:
	rm -f $(pspc_tests_OBJS) $(pspc_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	rm -f log count 
ifeq ($(SRC_DIR),$(BLD_DIR))
	-cd field; $(MAKE) clean
	-cd solvers; $(MAKE) clean
	-cd system; $(MAKE) clean
	-cd sweep; $(MAKE) clean
	rm -f */Test */*.o */*.d 
endif

-include $(pspc_tests_OBJS:.o=.d)
-include $(pspc_OBJS:.o=.d)
