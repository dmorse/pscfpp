pspg_tests_=pspg/tests/Test.cu

pspg_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_tests_:.cu=.o))
