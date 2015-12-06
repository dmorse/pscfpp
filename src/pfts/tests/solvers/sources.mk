pfts_tests_solvers_= pfts/tests/solvers/Test.cc

pfts_tests_solvers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pfts_tests_solvers_))
pfts_tests_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pfts_tests_solvers_:.cc=.o))

