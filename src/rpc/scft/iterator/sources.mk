rpc_scft_iterator_= \
  rpc/scft/iterator/AmIteratorBasis.cpp \
  rpc/scft/iterator/IteratorFactory.cpp

rpc_scft_iterator_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_scft_iterator_:.cpp=.o))

