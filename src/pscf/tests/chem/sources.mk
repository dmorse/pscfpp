pscf_tests_chem_=pscf/tests/chem/Test.cc

pscf_tests_chem_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pscf_tests_chem_))
pscf_tests_chem_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_chem_:.cc=.o))

