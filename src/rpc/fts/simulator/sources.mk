rpc_fts_simulator_= \
  rpc/fts/simulator/Simulator.cpp \
  rpc/fts/simulator/SimulatorFactory.cpp \
  rpc/fts/simulator/SimState.cpp 
  
rpc_fts_simulator_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_fts_simulator_:.cpp=.o))

