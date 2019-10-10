pspg_iterator_= \
  pspg/iterator/Iterator.cu 

  

pspg_iterator_SRCS=\
	  $(addprefix $(SRC_DIR)/, $(pspg_iterator_))
pspg_iterator_OBJS=\
	  $(addprefix $(BLD_DIR)/, $(pspg_iterator_:.cu=.o))

