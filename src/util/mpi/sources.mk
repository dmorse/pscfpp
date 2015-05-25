
util_mpi_=util/mpi/MpiFileIo.cpp \
    util/mpi/MpiTraits.cpp 

ifdef UTIL_MPI
util_mpi_+=\
    util/mpi/MpiLogger.cpp \
    util/mpi/MpiSendRecv.cpp \
    util/mpi/MpiStructBuilder.cpp
endif

util_mpi_SRCS=$(addprefix $(SRC_DIR)/, $(util_mpi_))
util_mpi_OBJS=$(addprefix $(BLD_DIR)/, $(util_mpi_:.cpp=.o))

