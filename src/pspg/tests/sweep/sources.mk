pspg_tests_sweep_=pspg/tests/sweep/Test.cu

pspg_tests_sweep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_tests_sweep_:.cu=.o))

