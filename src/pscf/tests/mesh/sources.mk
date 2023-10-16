pscf_tests_mesh_=pscf/tests/mesh/Test.cpp

pscf_tests_mesh_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_mesh_:.cpp=.o))

