fd1d_= \
  fd1d/TridiagonalSolver.cpp

fd1d_SRCS=\
     $(addprefix $(SRC_DIR)/, $(fd1d_))
fd1d_OBJS=\
     $(addprefix $(BLD_DIR)/, $(fd1d_:.cpp=.o))

$(fd1d_LIB): $(fd1d_OBJS)
	$(AR) rcs $(fd1d_LIB) $(fd1d_OBJS)

