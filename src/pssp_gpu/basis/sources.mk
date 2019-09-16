pssp_gpu_basis_= \
  pssp_gpu/basis/Basis.cu

pssp_gpu_basis_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_gpu_basis_))
pssp_gpu_basis_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_gpu_basis_:.cu=.o))

