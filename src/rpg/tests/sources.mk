rpg_tests_=rpg/tests/Test.cu

rpg_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_tests_:.cu=.ou))

rpg_tests_DEPS=\
     $(addprefix $(BLD_DIR)/, $(rpg_tests_:.cu=.du))
