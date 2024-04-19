rpg_simulate_compressor_= \
  rpg/simulate/compressor/CompressorFactory.cu \
  rpg/simulate/compressor/AmCompressor.cu \
  rpg/simulate/compressor/LrAmCompressor.cu

rpg_simulate_compressor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_simulate_compressor_:.cu=.o))

