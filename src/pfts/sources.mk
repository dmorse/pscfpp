pfts_= \
  pfts/Monomer.cpp\
  pfts/Block.cpp\
  pfts/Vertex.cpp\
  pfts/PolymerDescriptor.cpp\
  pfts/Species.cpp

pfts_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pfts_))
pfts_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pfts_:.cpp=.o))

$(pfts_LIB): $(pfts_OBJS)
	$(AR) rcs $(pfts_LIB) $(pfts_OBJS)

