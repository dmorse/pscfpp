pfts_chem_= \
  pfts/chem/Monomer.cpp\
  pfts/chem/Block.cpp\
  pfts/chem/Vertex.cpp

pfts_chem_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pfts_chem_))
pfts_chem_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pfts_chem_:.cpp=.o))

