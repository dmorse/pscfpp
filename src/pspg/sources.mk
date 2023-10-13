include $(SRC_DIR)/pspg/solvers/sources.mk
include $(SRC_DIR)/pspg/iterator/sources.mk
include $(SRC_DIR)/pspg/sweep/sources.mk
include $(SRC_DIR)/pspg/compressor/sources.mk
include $(SRC_DIR)/pspg/simulate/sources.mk

pspg_= \
  $(pspg_solvers_) \
  $(pspg_iterator_) \
  $(pspg_sweep_) \
  $(pspg_compressor_) \
  $(pspg_simulate_) \
  pspg/System.cu

pspg_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_))
pspg_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_:.cu=.o))

$(pspg_LIB): $(pspg_OBJS)
	$(AR) rcs $(pspg_LIB) $(pspg_OBJS)

