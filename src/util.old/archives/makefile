SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/archives/sources.mk

all: $(util_archives_OBJS)

clean:
	rm -f $(util_archives_OBJS) $(util_archives_OBJS:.o=.d)

clean-deps:
	rm -f $(util_archives_OBJS:.o=.d)

-include $(util_archives_OBJS:.o=.d)

