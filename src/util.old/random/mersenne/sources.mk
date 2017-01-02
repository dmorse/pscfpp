
util_random_mersenne_=util/random/mersenne/mtrand.cpp

util_random_mersenne_SRCS=\
    $(addprefix $(SRC_DIR)/, $(util_random_mersenne_))
util_random_mersenne_OBJS=\
    $(addprefix $(BLD_DIR)/, $(util_random_mersenne_:.cpp=.o))

