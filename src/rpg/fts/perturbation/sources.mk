rpg_fts_perturbation_= \
  rpg/fts/perturbation/Perturbation.cu \
  rpg/fts/perturbation/PerturbationFactory.cu \
  rpg/fts/perturbation/EinsteinCrystalPerturbation.cu 
  
rpg_fts_perturbation_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_fts_perturbation_:.cu=.o))

