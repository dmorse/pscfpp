cyln_misc_= \
  cyln/misc/Domain.cpp

cyln_misc_SRCS=\
     $(addprefix $(SRC_DIR)/, $(cyln_misc_))
cyln_misc_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cyln_misc_:.cpp=.o))

