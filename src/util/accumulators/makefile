SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/accumulators/sources.mk

all: $(util_accumulators_OBJS)

clean:
	rm -f $(util_accumulators_OBJS) $(util_accumulators_OBJS:.o=.d)

clean-deps:
	rm -f $(util_accumulators_OBJS:.o=.d)

-include $(util_accumulators_OBJS:.o=.d)

