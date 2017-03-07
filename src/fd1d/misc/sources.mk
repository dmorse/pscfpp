
fd1d_misc_=\
  fd1d/misc/HomogeneousComparison.cpp \
  fd1d/misc/FieldIo.cpp 

fd1d_misc_SRCS=\
     $(addprefix $(SRC_DIR)/, $(fd1d_misc_))
fd1d_misc_OBJS=\
     $(addprefix $(BLD_DIR)/, $(fd1d_misc_:.cpp=.o))

