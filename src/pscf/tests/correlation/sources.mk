pscf_tests_correlation_=pscf/tests/correlation/Test.cpp

pscf_tests_correlation_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_correlation_:.cpp=.o))

pscf_tests_correlation_DEPS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_correlation_:.cpp=.d))

