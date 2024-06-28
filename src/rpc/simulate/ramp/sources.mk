rpc_simulate_ramp_= \
  rpc/simulate/ramp/Ramp.cpp \
  rpc/simulate/ramp/LinearRamp.cpp \
  rpc/simulate/ramp/RampFactory.cpp 
  
rpc_simulate_ramp_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_simulate_ramp_:.cpp=.o))

