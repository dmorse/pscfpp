rpg_tests_fts_=rpg/tests/fts/Test.cu

rpg_tests_fts_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_tests_fts_:.cu=.o))

