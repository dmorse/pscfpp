rpg_simulate_ramp_= \
  rpg/simulate/ramp/Ramp.cu \
  rpg/simulate/ramp/LinearRamp.cu \
  rpg/simulate/ramp/RampFactory.cu 
  
rpg_simulate_ramp_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_simulate_ramp_:.cu=.o))

