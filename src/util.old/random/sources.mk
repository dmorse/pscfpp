include $(SRC_DIR)/util/random/mersenne/sources.mk

util_random_=$(util_random_mersenne_) \
    util/random/Ar1Process.cpp \
    util/random/Random.cpp 

util_random_SRCS=$(addprefix $(SRC_DIR)/, $(util_random_))
util_random_OBJS=$(addprefix $(BLD_DIR)/, $(util_random_:.cpp=.o))

