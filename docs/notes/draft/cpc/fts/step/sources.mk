cpc_fts_step_= \
  cpc/fts/step/Step.cpp \
  cpc/fts/step/StepFactory.cpp 
  
cpc_fts_step_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cpc_fts_step_:.cpp=.o))
