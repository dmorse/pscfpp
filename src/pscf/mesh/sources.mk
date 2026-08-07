pscf_mesh_CPP= \
  pscf/mesh/Mesh.cpp \
  pscf/mesh/MeshIterator.cpp \
  pscf/mesh/MeshIteratorFortran.cpp 

pscf_mesh_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_mesh_CPP:.cpp=.o))

pscf_mesh_DEPS=\
     $(addprefix $(BLD_DIR)/, $(pscf_mesh_CPP:.cpp=.d))

