include $(SRC_DIR)/rpg/field/sources.mk
include $(SRC_DIR)/rpg/solvers/sources.mk
include $(SRC_DIR)/rpg/iterator/sources.mk
include $(SRC_DIR)/rpg/sweep/sources.mk
include $(SRC_DIR)/rpg/compressor/sources.mk
include $(SRC_DIR)/rpg/simulate/sources.mk

rpg_= \
  $(rpg_field_) \
  $(rpg_solvers_) \
  $(rpg_iterator_) \
  $(rpg_sweep_) \
  $(rpg_compressor_) \
  $(rpg_simulate_) \
  rpg/System.cu

rpg_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_:.cu=.o))

$(rpg_LIB): $(rpg_OBJS)
	$(AR) rcs $(rpg_LIB) $(rpg_OBJS)

