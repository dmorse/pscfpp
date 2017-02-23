pscf_tests_crystal_=pscf/tests/crystal/Test.cc

pscf_tests_crystal_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pscf_tests_crystal_))
pscf_tests_crystal_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_crystal_:.cc=.o))

