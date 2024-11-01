rpg_fts_mcmove_= \
  rpg/fts/montecarlo/McSimulator.cu \
  rpg/fts/montecarlo/McMove.cu \
  rpg/fts/montecarlo/McMoveFactory.cu \
  rpg/fts/montecarlo/McMoveManager.cu \
  rpg/fts/montecarlo/RealMove.cu \
  rpg/fts/montecarlo/FourierMove.cu \
  rpg/fts/montecarlo/ForceBiasMove.cu 
  
rpg_fts_mcmove_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_fts_mcmove_:.cu=.o))

