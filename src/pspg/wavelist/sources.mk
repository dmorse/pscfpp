pspg_wavelist_= \
  pspg/wavelist/WaveList.cu

pspg_wavelist_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_wavelist_))
pspg_wavelist_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_wavelist_:.cu=.o))

