rpc_iterator_= \
  rpc/iterator/AmIteratorBasis.cpp \
  rpc/iterator/FilmIterator.cpp \
  rpc/iterator/ImposedFieldsGenerator.cpp \
  rpc/iterator/IteratorFactory.cpp \
  rpc/iterator/MaskGenFilm.cpp

rpc_iterator_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_iterator_:.cpp=.o))

