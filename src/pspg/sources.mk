include $(SRC_DIR)/pspg/wavelist/sources.mk
include $(SRC_DIR)/pspg/field/sources.mk
include $(SRC_DIR)/pspg/iterator/sources.mk
include $(SRC_DIR)/pspg/solvers/sources.mk
include $(SRC_DIR)/pspg/crystal/sources.mk

pspg_= \
  $(pspg_wavelist_) \
  $(pspg_field_) \
  $(pspg_solvers_) \
  $(pspg_iterator_)\
  $(pspg_crystal_)

pspg_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_))
pspg_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_:.cu=.o))

$(pspg_LIB): $(pspg_OBJS)
	$(AR) rcs $(pspg_LIB) $(pspg_OBJS)

