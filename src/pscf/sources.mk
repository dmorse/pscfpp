pscf_= \
  pscf/Monomer.cpp \
  pscf/Vertex.cpp \
  pscf/BlockDescriptor.cpp \
  pscf/Species.cpp \
  pscf/TridiagonalSolver.cpp

pscf_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pscf_))
pscf_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_:.cpp=.o))

$(pscf_LIB): $(pscf_OBJS)
	$(AR) rcs $(pscf_LIB) $(pscf_OBJS)

