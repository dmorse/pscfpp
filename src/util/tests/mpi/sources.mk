util_tests_mpi_=util/tests/mpi/Test.cc

util_tests_mpi_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_mpi_))
util_tests_mpi_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_mpi_:.cc=.o))

