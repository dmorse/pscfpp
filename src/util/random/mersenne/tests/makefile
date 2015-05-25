SRC_DIR_REL =../../../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include sources.mk

all: $(util_random_mersenne_tests_OBJS) mttest

clean:
	rm -f $(util_random_mersenne_tests_OBJS) $(util_random_mersenne_tests_OBJS:.o=.d) \
        mttest

clean-deps:
	rm -f $(util_random_mersenne_tests_OBJS:.o=.d)

mttest: mttest.o
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o mttest mttest.o  ../mtrand.o

-include $(util_random_mersenne_tests_OBJS:.o=.d)

