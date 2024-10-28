rpg_fts_ramp_= \
  rpg/fts/ramp/Ramp.cu \
  rpg/fts/ramp/LinearRamp.cu \
  rpg/fts/ramp/RampFactory.cu 
  
rpg_fts_ramp_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_fts_ramp_:.cu=.o))

