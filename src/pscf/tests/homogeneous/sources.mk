pscf_tests_homogeneous_=pscf/tests/homogeneous/Test.cc

pscf_tests_homogeneous_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pscf_tests_homogeneous_))
pscf_tests_homogeneous_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_homogeneous_:.cc=.o))

