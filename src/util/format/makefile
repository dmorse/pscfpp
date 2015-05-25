SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/format/sources.mk

all: $(util_format_OBJS)

clean:
	rm -f $(util_format_OBJS) $(util_format_OBJS:.o=.d)

clean-deps:
	rm -f $(util_format_OBJS:.o=.d)

-include $(util_format_OBJS:.o=.d)

