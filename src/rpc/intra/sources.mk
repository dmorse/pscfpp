rpc_intra_= \
  rpc/intra/IntraCorrelation.cpp 

rpc_intra_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_intra_:.cpp=.o))

