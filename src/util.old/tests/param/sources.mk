include $(SRC_DIR)/util/param/tests/serial/sources.mk
include $(SRC_DIR)/util/param/tests/mpi/sources.mk

tests_param_util_SRCS=$(tests_param_util_serial_SRCS) \
    $(tests_param_util_mpi_SRCS) 

tests_param_util_OBJS=$(tests_param_util_SRCS:.cc=.o)

