pssp_gpu_iterator_= \
  pssp_gpu/iterator/Iterator.cu 

  

pssp_gpu_iterator_SRCS=\
	  $(addprefix $(SRC_DIR)/, $(pssp_gpu_iterator_))
pssp_gpu_iterator_OBJS=\
	  $(addprefix $(BLD_DIR)/, $(pssp_gpu_iterator_:.cu=.o))

