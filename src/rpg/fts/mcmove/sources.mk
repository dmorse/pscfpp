rpg_fts_mcmove_= \
  rpg/fts/mcmove/McSimulator.cu \
  rpg/fts/mcmove/McMove.cu \
  rpg/fts/mcmove/McMoveFactory.cu \
  rpg/fts/mcmove/McMoveManager.cu \
  rpg/fts/mcmove/RealMove.cu \
  rpg/fts/mcmove/FourierMove.cu \
  rpg/fts/mcmove/ForceBiasMove.cu 
  
rpg_fts_mcmove_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_fts_mcmove_:.cu=.o))

