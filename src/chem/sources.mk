chem_= \
  chem/Monomer.cpp\
  chem/Block.cpp\
  chem/Vertex.cpp

chem_SRCS=\
     $(addprefix $(SRC_DIR)/, $(chem_))
chem_OBJS=\
     $(addprefix $(BLD_DIR)/, $(chem_:.cpp=.o))

$(chem_LIB): $(chem_OBJS)
	$(AR) rcs $(chem_LIB) $(chem_OBJS)

