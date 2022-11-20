pspg_iterator_= \
  pspg/iterator/IteratorFactory.cu \
  pspg/iterator/AmIteratorBasis.cu \
  pspg/iterator/AmIteratorGrid.cu

pspg_iterator_SRCS=\
	  $(addprefix $(SRC_DIR)/, $(pspg_iterator_))
pspg_iterator_OBJS=\
	  $(addprefix $(BLD_DIR)/, $(pspg_iterator_:.cu=.o))

