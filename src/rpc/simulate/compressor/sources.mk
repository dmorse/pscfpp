# Include source files lists from subdirectories
include $(SRC_DIR)/rpc/simulate/compressor/intra/sources.mk

rpc_simulate_compressor_= \
  $(rpc_simulate_compressor_intra_) \
  rpc/simulate/compressor/CompressorFactory.cpp \
  rpc/simulate/compressor/AmCompressor.cpp \
  rpc/simulate/compressor/LrAmCompressor.cpp \
  rpc/simulate/compressor/LrCompressor.cpp \
  rpc/simulate/compressor/LrPostAmCompressor.cpp 

rpc_simulate_compressor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_simulate_compressor_:.cpp=.o))

