pscf_tests_homogeneous_=pscf/tests/homogeneous/Test.cpp

pscf_tests_homogeneous_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_homogeneous_:.cpp=.o))

