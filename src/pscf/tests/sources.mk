pscf_tests_=pscf/tests/Test.cpp
pscf_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_:.cpp=.o))

