rpg_fts_montecarlo_= \
  rpg/fts/montecarlo/McSimulator.cu \
  rpg/fts/montecarlo/McMove.cu \
  rpg/fts/montecarlo/McMoveFactory.cu \
  rpg/fts/montecarlo/McMoveManager.cu \
  rpg/fts/montecarlo/RealMove.cu \
  rpg/fts/montecarlo/ForceBiasMove.cu \
  rpg/fts/montecarlo/BdMove.cu \
  rpg/fts/montecarlo/ShiftMove.cu
  
rpg_fts_montecarlo_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_fts_montecarlo_:.cu=.ou))

