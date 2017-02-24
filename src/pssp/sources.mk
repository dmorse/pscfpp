pssp_= 
  # pssp/Basis.cpp 

pssp_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_))
pssp_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_:.cpp=.o))

$(pssp_LIB): $(pssp_OBJS)
	$(AR) rcs $(pssp_LIB) $(pssp_OBJS)

