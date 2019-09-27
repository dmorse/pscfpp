pssp_gpu_crystal_= \
  pssp_gpu/crystal/UnitCell1.cu \
  pssp_gpu/crystal/UnitCell2.cu \
  pssp_gpu/crystal/UnitCell3.cu \
  pssp_gpu/crystal/shiftToMinimum.cu

pssp_gpu_crystal_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_gpu_crystal_))
pssp_gpu_crystal_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_gpu_crystal_:.cu=.o))

