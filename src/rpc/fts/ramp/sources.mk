rpc_fts_ramp_= \
  rpc/fts/ramp/Ramp.cpp \
  rpc/fts/ramp/LinearRamp.cpp \
  rpc/fts/ramp/RampFactory.cpp 
  
rpc_fts_ramp_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_fts_ramp_:.cpp=.o))

