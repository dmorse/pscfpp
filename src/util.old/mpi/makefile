SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/mpi/sources.mk

all: $(util_mpi_OBJS)

clean:
	rm -f $(util_mpi_OBJS) $(util_mpi_OBJS:.o=.d)

clean-deps:
	rm -f $(util_mpi_OBJS:.o=.d)

-include $(util_mpi_OBJS:.o=.d)

