fd1d_= \
  fd1d/Propagator.cpp \
  fd1d/Block.cpp \
  fd1d/Polymer.cpp \
  fd1d/Solvent.cpp \
  fd1d/Mixture.cpp \
  fd1d/Iterator.cpp \
  fd1d/NrIterator.cpp 

fd1d_SRCS=\
     $(addprefix $(SRC_DIR)/, $(fd1d_))
fd1d_OBJS=\
     $(addprefix $(BLD_DIR)/, $(fd1d_:.cpp=.o))

$(fd1d_LIB): $(fd1d_OBJS)
	$(AR) rcs $(fd1d_LIB) $(fd1d_OBJS)

