util_math_=util/math/Constants.cpp 

util_math_SRCS=$(addprefix $(SRC_DIR)/, $(util_math_))
util_math_OBJS=$(addprefix $(BLD_DIR)/, $(util_math_:.cpp=.o))

