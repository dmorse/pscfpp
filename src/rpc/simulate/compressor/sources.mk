rpc_simulate_compressor_= \
  rpc/simulate/compressor/CompressorFactory.cpp \
  rpc/simulate/compressor/AmCompressor.cpp \
  rpc/simulate/compressor/LrAmCompressor.cpp \
  rpc/simulate/compressor/LrCompressor.cpp

rpc_simulate_compressor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_simulate_compressor_:.cpp=.o))

