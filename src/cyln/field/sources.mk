cyln_field_= \
  cyln/field/FFT.cpp 

cyln_field_SRCS=\
     $(addprefix $(SRC_DIR)/, $(cyln_field_))
cyln_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cyln_field_:.cpp=.o))

