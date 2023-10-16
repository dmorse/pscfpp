pspc_tests_sweep_=pspc/tests/sweep/Test.cpp

pspc_tests_sweep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_tests_sweep_:.cpp=.o))

