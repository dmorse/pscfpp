rpg_simulate_mcmove_= \
  rpg/simulate/mcmove/McSimulator.cu \
  rpg/simulate/mcmove/McState.cu \
  rpg/simulate/mcmove/McMove.cu \
  rpg/simulate/mcmove/McMoveFactory.cu \
  rpg/simulate/mcmove/McMoveManager.cu \
  rpg/simulate/mcmove/RealMove.cu \
  rpg/simulate/mcmove/FourierMove.cu \
  rpg/simulate/mcmove/ForceBiasMove.cu 
  
rpg_simulate_mcmove_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_simulate_mcmove_:.cu=.o))

