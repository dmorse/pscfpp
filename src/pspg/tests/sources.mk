pspg_tests_=pspg/tests/Test.ccu

pspg_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_tests_))
pspg_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_tests_:.cu=.o))
