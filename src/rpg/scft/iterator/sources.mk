rpg_scft_iterator_= \
  rpg/scft/iterator/IteratorFactory.cu \
  rpg/scft/iterator/AmIteratorBasis.cu \
  rpg/scft/iterator/AmIteratorGrid.cu

rpg_scft_iterator_OBJS=\
	  $(addprefix $(BLD_DIR)/, $(rpg_scft_iterator_:.cu=.o))

