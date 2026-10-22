r1_iterator_CPP=\
  r1/iterator/Iterator.cpp \
  r1/iterator/NrIterator.cpp \
  r1/iterator/BinaryRelaxIterator.cpp \
  r1/iterator/AmIterator.cpp \
  r1/iterator/IteratorFactory.cpp

r1_iterator_OBJS=\
     $(addprefix $(BLD_DIR)/, $(r1_iterator_CPP:.cpp=.o))

r1_iterator_DEPS=\
     $(addprefix $(BLD_DIR)/, $(r1_iterator_CPP:.cpp=.d))
