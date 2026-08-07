pscf_tests_chem_=pscf/tests/chem/Test.cpp

pscf_tests_chem_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_chem_:.cpp=.o))

pscf_tests_chem_DEPS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_chem_:.cpp=.d))

