cpc_tests_solvers_=cpc/tests/solvers/Test.cpp

cpc_tests_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cpc_tests_solvers_:.cpp=.o))

