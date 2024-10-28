rpg_fts_compressor_intra_= \
  rpg/fts/compressor/intra/IntraCorrelation.cu 

rpg_fts_compressor_intra_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_fts_compressor_intra_:.cu=.o))

