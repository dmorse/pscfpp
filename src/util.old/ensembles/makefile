SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/ensembles/sources.mk

all: $(util_ensembles_OBJS)

clean:
	rm -f $(util_ensembles_OBJS) $(util_ensembles_OBJS:.o=.d)

clean-deps:
	rm -f $(util_ensembles_OBJS:.o=.d)

-include $(util_ensembles_OBJS:.o=.d)

