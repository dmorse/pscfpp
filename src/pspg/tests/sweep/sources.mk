pspg_tests_sweep_=pspg/tests/sweep/Test.cc

pspg_tests_sweep_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_tests_sweep_))
pspg_tests_sweep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_tests_sweep_:.cc=.o))

