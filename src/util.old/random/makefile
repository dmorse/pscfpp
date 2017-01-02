SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/random/sources.mk

all: $(util_random_OBJS)

clean:
	rm -f $(util_random_OBJS) $(util_random_OBJS:.o=.d)

clean-deps:
	rm -f $(util_random_OBJS:.o=.d)

-include $(util_random_OBJS:.o=.d)

