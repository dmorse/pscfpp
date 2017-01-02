
util_boundary_=\
    util/boundary/OrthoRegion.cpp \
    util/boundary/OrthorhombicBoundary.cpp \
    util/boundary/MonoclinicBoundary.cpp 
    #util/boundary/MonoclinicBoundaryMI.cpp

util_boundary_SRCS=$(addprefix $(SRC_DIR)/, $(util_boundary_))
util_boundary_OBJS=$(addprefix $(BLD_DIR)/, $(util_boundary_:.cpp=.o))

