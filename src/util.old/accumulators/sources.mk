
util_accumulators_=util/accumulators/Average.cpp \
    util/accumulators/AverageStage.cpp \
    util/accumulators/TensorAverage.cpp \
    util/accumulators/SymmTensorAverage.cpp \
    util/accumulators/Distribution.cpp \
    util/accumulators/IntDistribution.cpp \
    util/accumulators/RadialDistribution.cpp 

util_accumulators_SRCS=$(addprefix $(SRC_DIR)/, $(util_accumulators_))
util_accumulators_OBJS=$(addprefix $(BLD_DIR)/, $(util_accumulators_:.cpp=.o))

