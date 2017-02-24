pscf_tests_mesh_=pscf/tests/mesh/Test.cc

pscf_tests_mesh_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pscf_tests_mesh_))
pscf_tests_mesh_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_mesh_:.cc=.o))

