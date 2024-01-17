pspg_tests_simulate_=pspg/tests/simulate/Test.cu

pspg_tests_simulate_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_tests_simulate_:.cu=.o))

