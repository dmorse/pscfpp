
fd1d_iterator_=\
  fd1d/iterator/Iterator.cpp 

ifdef PSCF_GSL
fd1d_iterator_+= fd1d/iterator/NrIterator.cpp
endif

fd1d_iterator_SRCS=\
     $(addprefix $(SRC_DIR)/, $(fd1d_iterator_))
fd1d_iterator_OBJS=\
     $(addprefix $(BLD_DIR)/, $(fd1d_iterator_:.cpp=.o))

