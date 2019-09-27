include $(SRC_DIR)/pssp_gpu/basis/sources.mk
include $(SRC_DIR)/pssp_gpu/field/sources.mk
include $(SRC_DIR)/pssp_gpu/iterator/sources.mk
include $(SRC_DIR)/pssp_gpu/solvers/sources.mk
include $(SRC_DIR)/pssp_gpu/crystal/sources.mk
include $(SRC_DIR)/pssp_gpu/inter/sources.mk

pssp_gpu_= \
  $(pssp_gpu_basis_) \
  $(pssp_gpu_field_) \
  $(pssp_gpu_solvers_) \
  $(pssp_gpu_iterator_)\
  $(pssp_gpu_crystal_)\
  $(pssp_gpu_inter_)


pssp_gpu_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_gpu_))
pssp_gpu_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_gpu_:.cu=.o))

$(pssp_gpu_LIB): $(pssp_gpu_OBJS)
	$(AR) rcs $(pssp_gpu_LIB) $(pssp_gpu_OBJS)

