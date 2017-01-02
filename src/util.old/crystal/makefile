SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/crystal/sources.mk

all: $(util_crystal_OBJS)

clean:
	rm -f $(util_crystal_OBJS) $(util_crystal_OBJS:.o=.d)

clean-deps:
	rm -f $(util_crystal_OBJS:.o=.d)

-include $(util_crystal_OBJS:.o=.d)

