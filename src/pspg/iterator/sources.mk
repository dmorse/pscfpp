pspg_iterator_= \
  pspg/iterator/AmIterator.cu 

  

pspg_iterator_SRCS=\
	  $(addprefix $(SRC_DIR)/, $(pspg_iterator_))
pspg_iterator_OBJS=\
	  $(addprefix $(BLD_DIR)/, $(pspg_iterator_:.cu=.o))

