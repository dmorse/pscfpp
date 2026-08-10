# -------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/rp/sources.mk
include $(SRC_DIR)/rp/patterns.mk

# -------------------------------------------------------------------------
# Main targets

all: Test.o Test

Test: Test.o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o Test Test.o $(LIBS)

run: Test.o Test
	./Test $(SRC_DIR)/rp/tests/cpp/ 

run-quiet: Test.o Test
	@./Test $(SRC_DIR)/rp/tests/cpp/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in rp/tests/cpp/log" > count
	@cat count

clean-outputs:
	rm -f log count

clean:
	rm -f Test
	rm -f Test.o Test.d
	$(MAKE) clean-outputs

.PHONY= all run run-quiet clean-outputs clean
# -------------------------------------------------------------------------
# Main targets

-include Test.d
-include $(rp_DEPS)
-include $(prdc_DEPS)
-include $(pscf_DEPS)
-include $(util_DEPS)
