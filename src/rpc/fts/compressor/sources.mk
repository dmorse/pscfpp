# Include source files lists from subdirectories
include $(SRC_DIR)/rpc/fts/compressor/intra/sources.mk

rpc_fts_compressor_= \
  $(rpc_fts_compressor_intra_) \
  rpc/fts/compressor/CompressorFactory.cpp \
  rpc/fts/compressor/AmCompressor.cpp \
  rpc/fts/compressor/LrCompressor.cpp \
  rpc/fts/compressor/LrAmCompressor.cpp 

rpc_fts_compressor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_fts_compressor_:.cpp=.o))

