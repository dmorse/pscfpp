rpg_tests_environment_=rpg/tests/environment/Test.cu

rpg_tests_environment_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_tests_environment_:.cu=.o))

