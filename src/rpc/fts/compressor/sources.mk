rpc_fts_compressor_= \
  rpc/fts/compressor/CompressorFactory.cpp \
  rpc/fts/compressor/AmCompressor.cpp \
  rpc/fts/compressor/LrCompressor.cpp \
  rpc/fts/compressor/LrAmCompressor.cpp \
  rpc/fts/compressor/IntraCorrelation.cpp

rpc_fts_compressor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_fts_compressor_:.cpp=.o))

