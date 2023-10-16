pscf_tests_inter_=pscf/tests/inter/Test.cpp

pscf_tests_inter_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_inter_:.cpp=.o))

