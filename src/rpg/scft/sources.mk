# Include source files lists from subdirectories
include $(SRC_DIR)/rpg/scft/iterator/sources.mk
include $(SRC_DIR)/rpg/scft/sweep/sources.mk

rpg_scft_= \
  $(rpg_scft_iterator_) \
  $(rpg_scft_sweep_) \
  
rpg_scft_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_scft_:.cu=.o))

