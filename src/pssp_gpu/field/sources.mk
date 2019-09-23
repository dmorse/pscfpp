pssp_gpu_field_= \
  pssp_gpu/field/FFT.cu \
  pssp_gpu/field/FFTBatched.cu

pssp_gpu_field_SRCS=\
	  $(addprefix $(SRC_DIR)/, $(pssp_gpu_field_))
pssp_gpu_field_OBJS=\
	  $(addprefix $(BLD_DIR)/, $(pssp_gpu_field_:.cu=.o))

