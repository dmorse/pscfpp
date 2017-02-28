pssp_tests_solvers_=pssp/tests/solvers/Test.cc

pssp_tests_solvers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_tests_solvers_))
pssp_tests_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_tests_solvers_:.cc=.o))

