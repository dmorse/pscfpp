SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/containers/sources.mk

all: $(util_containers_OBJS)

clean:
	rm -f $(util_containers_OBJS) $(util_containers_OBJS:.o=.d)

clean-deps:
	rm -f $(util_containers_OBJS:.o=.d)

-include $(util_containers_OBJS:.o=.d)

