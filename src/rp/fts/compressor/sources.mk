
rp_fts_compressor_OBJS=
rp_fts_compressor_DEPS=

ifdef PSCF_CPP
  rp_fts_compressor_CPP= \
    rp/fts/compressor/Compressor.cpp \
    rp/fts/compressor/IntraCorrelation.cpp \
    rp/fts/compressor/AmCompressor.cpp \
    rp/fts/compressor/LrCompressor.cpp \
    rp/fts/compressor/LrAmCompressor.cpp \
    rp/fts/compressor/CompressorFactory.cpp
  rp_fts_compressor_OBJS+=\
       $(addprefix $(BLD_DIR)/, $(rp_fts_compressor_CPP:.cpp=.o))
  rp_fts_compressor_DEPS+=\
       $(addprefix $(BLD_DIR)/, $(rp_fts_compressor_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_fts_compressor_CU= \
    rp/fts/compressor/Compressor.cu \
    rp/fts/compressor/IntraCorrelation.cu \
    rp/fts/compressor/AmCompressor.cu \
    rp/fts/compressor/LrCompressor.cu \
    rp/fts/compressor/LrAmCompressor.cu \
    rp/fts/compressor/CompressorFactory.cu
  rp_fts_compressor_OBJS+=\
      $(addprefix $(BLD_DIR)/, $(rp_fts_compressor_CU:.cu=.ou))
  rp_fts_compressor_DEPS+=\
      $(addprefix $(BLD_DIR)/, $(rp_fts_compressor_CU:.cu=.du))
endif
