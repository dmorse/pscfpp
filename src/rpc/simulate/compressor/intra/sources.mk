rpc_simulate_compressor_intra_= \
  rpc/simulate/compressor/intra/IntraCorrelation.cpp 

rpc_simulate_compressor_intra_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_simulate_compressor_intra_:.cpp=.o))

