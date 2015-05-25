SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/signal/sources.mk

all: $(util_signal_OBJS)

clean:
	rm -f $(util_signal_OBJS) $(util_signal_OBJS:.o=.d)

clean-deps:
	rm -f $(util_signal_OBJS:.o=.d)

-include $(util_signal_OBJS:.o=.d)

