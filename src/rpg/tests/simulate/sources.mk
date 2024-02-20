rpg_tests_simulate_=rpg/tests/simulate/Test.cu

rpg_tests_simulate_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_tests_simulate_:.cu=.o))

