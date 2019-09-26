pssp_gpu_wavelist_= \
  pssp_gpu/wavelist/WaveList.cu

pssp_gpu_wavelist_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_gpu_wavelist_))
pssp_gpu_wavelist_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_gpu_wavelist_:.cu=.o))

