pspc_compressor_= \
  pspc/compressor/CompressorFactory.cpp \
  pspc/compressor/AmCompressor.cpp \
  pspc/compressor/LrAmCompressor.cpp \
  pspc/compressor/LrCompressor.cpp

pspc_compressor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_compressor_:.cpp=.o))

