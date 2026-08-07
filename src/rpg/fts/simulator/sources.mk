rpg_fts_simulator_=\
  rpg/fts/simulator/SimState.cu \
  rpg/fts/simulator/Simulator.cu \
  rpg/fts/simulator/SimulatorFactory.cu 
  
rpg_fts_simulator_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_fts_simulator_:.cu=.ou))

