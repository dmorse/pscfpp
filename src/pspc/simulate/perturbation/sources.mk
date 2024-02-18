pspc_simulate_perturbation_= \
  pspc/simulate/perturbation/Perturbation.cpp \
  pspc/simulate/perturbation/PerturbationFactory.cpp 
  
pspc_simulate_perturbation_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_simulate_perturbation_:.cpp=.o))

