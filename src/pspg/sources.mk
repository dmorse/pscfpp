include $(SRC_DIR)/pspg/field/sources.mk
include $(SRC_DIR)/pspg/iterator/sources.mk
include $(SRC_DIR)/pspg/solvers/sources.mk

pspg_= \
  $(pspg_field_) \
  $(pspg_solvers_) \
  $(pspg_iterator_) \
  pspg/System.cu

pspg_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_))
pspg_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_:.cu=.o))

$(pspg_LIB): $(pspg_OBJS)
	$(AR) rcs $(pspg_LIB) $(pspg_OBJS)

