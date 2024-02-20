rpg_iterator_= \
  rpg/iterator/IteratorFactory.cu \
  rpg/iterator/AmIteratorBasis.cu \
  rpg/iterator/AmIteratorGrid.cu

rpg_iterator_OBJS=\
	  $(addprefix $(BLD_DIR)/, $(rpg_iterator_:.cu=.o))

