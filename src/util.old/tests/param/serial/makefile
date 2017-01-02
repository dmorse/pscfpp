BLD_DIR_REL =../../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/util/tests/param/serial/sources.mk

all: $(util_tests_param_serial_OBJS)

clean:
	rm -f $(util_tests_param_serial_OBJS) 
	rm -f $(util_tests_param_serial_OBJS:.o=.d)
	rm -f $(util_tests_param_serial_OBJS:.o=)
	$(MAKE) clean-outputs

clean-outputs:
	rm -f out/*

-include $(util_tests_param_serial_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

