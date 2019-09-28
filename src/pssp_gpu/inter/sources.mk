pssp_gpu_inter_= \
  pssp_gpu/inter/ChiInteraction.cu 

pssp_gpu_inter_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_gpu_inter_))
pssp_gpu_inter_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_gpu_inter_:.cu=.o))

