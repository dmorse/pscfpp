rpg_simulate_perturbation_= \
  rpg/simulate/perturbation/Perturbation.cu \
  rpg/simulate/perturbation/PerturbationFactory.cu \
  rpg/simulate/perturbation/EinsteinCrystalPerturbation.cu 
  
rpg_simulate_perturbation_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_simulate_perturbation_:.cu=.o))

