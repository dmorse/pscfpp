BLD_DIR_REL =../..

include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/misc/sources.mk

all: $(util_misc_OBJS)

clean:
	rm -f $(util_misc_OBJS) $(util_misc_OBJS:.o=.d)

clean-deps:
	rm -f $(util_misc_OBJS:.o=.d)

-include $(util_misc_OBJS:.o=.d)

