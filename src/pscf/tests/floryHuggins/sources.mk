pscf_tests_floryHuggins_=pscf/tests/floryHuggins/Test.cpp

pscf_tests_floryHuggins_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_floryHuggins_:.cpp=.o))

pscf_tests_floryHuggins_DEPS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_floryHuggins_:.cpp=.d))

