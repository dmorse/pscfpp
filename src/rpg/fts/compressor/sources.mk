# Include source files lists from subdirectories
include $(SRC_DIR)/rpg/fts/compressor/intra/sources.mk

rpg_fts_compressor_= \
  $(rpg_fts_compressor_intra_) \
  rpg/fts/compressor/CompressorFactory.cu \
  rpg/fts/compressor/AmCompressor.cu \
  rpg/fts/compressor/LrCompressor.cu  \
  rpg/fts/compressor/LrPostAmCompressor.cu \
  rpg/fts/compressor/LrAmPreCompressor.cu 
  

rpg_fts_compressor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_fts_compressor_:.cu=.o))

