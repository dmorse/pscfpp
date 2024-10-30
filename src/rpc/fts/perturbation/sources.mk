rpc_fts_perturbation_= \
  rpc/fts/perturbation/Perturbation.cpp \
  rpc/fts/perturbation/PerturbationFactory.cpp \
  rpc/fts/perturbation/EinsteinCrystalPerturbation.cpp 
  
rpc_fts_perturbation_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_fts_perturbation_:.cpp=.o))

