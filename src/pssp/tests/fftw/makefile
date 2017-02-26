BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/pssp/include.mk
include $(SRC_DIR)/pssp/tests/fftw/sources.mk

TEST=pssp/tests/fftw/Test

all: $(pssp_tests_fftw_OBJS) $(BLD_DIR)/$(TEST)

includes:
	@echo $(INCLUDES)

libs:
	@echo $(LIBS)

run: $(pssp_tests_fftw_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pssp/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pssp/tests/log" > count
	@cat count

clean:
	rm -f $(pssp_tests_fftw_OBJS) $(pssp_tests_fftw_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	rm -f log count 

-include $(pssp_tests_fftw_OBJS:.o=.d)
-include $(pssp_tests_fftw_OBJS:.o=.d)
