
fd1d_commands_=\
  fd1d/commands/HomogeneousCommand.cpp 

fd1d_commands_SRCS=\
     $(addprefix $(SRC_DIR)/, $(fd1d_commands_))
fd1d_commands_OBJS=\
     $(addprefix $(BLD_DIR)/, $(fd1d_commands_:.cpp=.o))

