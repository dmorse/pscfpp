SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/boundary/sources.mk

all: $(util_boundary_OBJS)

clean:
	rm -f $(util_boundary_OBJS) $(util_boundary_OBJS:.o=.d)

clean-deps:
	rm -f $(util_boundary_OBJS:.o=.d)

-include $(util_boundary_OBJS:.o=.d)

