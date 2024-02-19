rpc_iterator_= \
  rpc/iterator/AmIteratorBasis.cpp \
  rpc/iterator/FilmIterator.cpp \
  rpc/iterator/IteratorFactory.cpp 

rpc_iterator_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_iterator_:.cpp=.o))

