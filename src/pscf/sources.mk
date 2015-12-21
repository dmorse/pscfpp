pscf_= \
  pscf/Monomer.cpp \
  pscf/Block.cpp \
  pscf/Vertex.cpp \
  pscf/PolymerDescriptor.cpp \
  pscf/Species.cpp \
  pscf/TridiagonalSolver.cpp

pscf_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pscf_))
pscf_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_:.cpp=.o))

$(pscf_LIB): $(pscf_OBJS)
	$(AR) rcs $(pscf_LIB) $(pscf_OBJS)

