pspc_tests_simulate_=pspc/tests/simulate/Test.cpp

pspc_tests_simulate_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_tests_simulate_:.cpp=.o))

