pspg_tests_solvers_=pspg/tests/solvers/Test.ccu

pspg_tests_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_tests_solvers_:.ccu=.o))

