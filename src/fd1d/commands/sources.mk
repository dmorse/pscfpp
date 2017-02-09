
fd1d_commands_=\
  fd1d/commands/HomogeneousComparison.cpp \
  fd1d/commands/FieldEditor.cpp 

fd1d_commands_SRCS=\
     $(addprefix $(SRC_DIR)/, $(fd1d_commands_))
fd1d_commands_OBJS=\
     $(addprefix $(BLD_DIR)/, $(fd1d_commands_:.cpp=.o))

