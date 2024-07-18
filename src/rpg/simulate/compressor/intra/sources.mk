rpg_simulate_compressor_intra_= \
  rpg/simulate/compressor/intra/IntraCorrelation.cu 

rpg_simulate_compressor_intra_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_simulate_compressor_intra_:.cu=.o))

