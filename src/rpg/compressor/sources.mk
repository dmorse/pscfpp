rpg_compressor_= \
  rpg/compressor/CompressorFactory.cu \
  rpg/compressor/AmCompressor.cu \
  rpg/compressor/LrAmCompressor.cu

rpg_compressor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_compressor_:.cu=.o))

