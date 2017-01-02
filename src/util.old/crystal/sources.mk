
util_crystal_=\
    util/crystal/LatticeSystem.cpp \
    util/crystal/PointGroup.cpp \
    util/crystal/PointSymmetry.cpp 

util_crystal_SRCS=$(addprefix $(SRC_DIR)/, $(util_crystal_))
util_crystal_OBJS=$(addprefix $(BLD_DIR)/, $(util_crystal_:.cpp=.o))

