util_tests_param_mpi_=util/tests/param/mpi/MpiTest.cc


util_tests_param_mpi_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_param_mpi_))
util_tests_param_mpi_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_param_mpi_:.cc=.o))

