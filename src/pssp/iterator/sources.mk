pssp_iterator_= \
  pssp/iterator/Iterator.cpp \
  pssp/iterator/AmIterator.cpp

  

pssp_iterator_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_iterator_))
pssp_iterator_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_iterator_:.cpp=.o))

