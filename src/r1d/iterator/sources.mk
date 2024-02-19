r1d_iterator_=\
  r1d/iterator/Iterator.cpp \
  r1d/iterator/NrIterator.cpp \
  r1d/iterator/BinaryRelaxIterator.cpp \
  r1d/iterator/AmIterator.cpp \
  r1d/iterator/IteratorFactory.cpp

r1d_iterator_OBJS=\
     $(addprefix $(BLD_DIR)/, $(r1d_iterator_:.cpp=.o))
