pspg_tests_solvers_=pspg/tests/solvers/Test.cu

pspg_tests_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_tests_solvers_:.cu=.o))

