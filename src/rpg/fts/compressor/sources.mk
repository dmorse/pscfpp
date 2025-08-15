rpg_fts_compressor_= \
  rpg/fts/compressor/CompressorFactory.cu \
  rpg/fts/compressor/AmCompressor.cu \
  rpg/fts/compressor/LrCompressor.cu  \
  rpg/fts/compressor/LrAmCompressor.cu \
  rpg/fts/compressor/IntraCorrelation.cu 

rpg_fts_compressor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_fts_compressor_:.cu=.o))

