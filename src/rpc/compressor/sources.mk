rpc_compressor_= \
  rpc/compressor/CompressorFactory.cpp \
  rpc/compressor/AmCompressor.cpp \
  rpc/compressor/LrAmCompressor.cpp \
  rpc/compressor/LrCompressor.cpp

rpc_compressor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_compressor_:.cpp=.o))

