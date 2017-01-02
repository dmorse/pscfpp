SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/param/sources.mk

all: $(util_param_OBJS)

clean:
	rm -f $(util_param_OBJS) $(util_param_OBJS:.o=.d)

clean-deps:
	rm -f $(util_param_OBJS:.o=.d)

-include $(util_param_OBJS:.o=.d)

