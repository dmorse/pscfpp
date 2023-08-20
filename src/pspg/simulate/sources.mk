pspg_simulate_= \
  pspg/simulate/McState.cu 
  
  
pspg_simulate_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_simulate_))
pspg_simulate_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_simulate_:.cu=.o))

