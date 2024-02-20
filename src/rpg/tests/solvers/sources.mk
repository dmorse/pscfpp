rpg_tests_solvers_=rpg/tests/solvers/Test.cu

rpg_tests_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_tests_solvers_:.cu=.o))

