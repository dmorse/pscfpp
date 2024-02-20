rpg_tests_system_=rpg/tests/system/Test.cu

rpg_tests_system_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_tests_system_:.cu=.o))

