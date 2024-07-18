# Include source files lists from subdirectories
include $(SRC_DIR)/rpg/simulate/compressor/intra/sources.mk

rpg_simulate_compressor_= \
  $(rpg_simulate_compressor_intra_) \
  rpg/simulate/compressor/CompressorFactory.cu \
  rpg/simulate/compressor/AmCompressor.cu \
  rpg/simulate/compressor/LrAmCompressor.cu \
  rpg/simulate/compressor/LrCompressor.cu  \
  rpg/simulate/compressor/LrPostAmCompressor.cu \
  

rpg_simulate_compressor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_simulate_compressor_:.cu=.o))

