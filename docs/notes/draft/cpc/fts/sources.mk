include $(SRC_DIR)/cpc/fts/step/sources.mk

cpc_fts_= \
  $(cpc_fts_step_) \
  cpc/fts/Simulator.cpp
  
cpc_fts_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cpc_fts_:.cpp=.o))
