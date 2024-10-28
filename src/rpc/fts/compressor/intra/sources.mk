rpc_fts_compressor_intra_= \
  rpc/fts/compressor/intra/IntraCorrelation.cpp 

rpc_fts_compressor_intra_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_fts_compressor_intra_:.cpp=.o))

