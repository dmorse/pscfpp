include $(SRC_DIR)/cyln/field/sources.mk
include $(SRC_DIR)/cyln/solvers/sources.mk

cyln_=\
  $(cyln_field_) \
  $(cyln_solvers_) 

cyln_SRCS=\
     $(addprefix $(SRC_DIR)/, $(cyln_))
cyln_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cyln_:.cpp=.o))

$(cyln_LIB): $(cyln_OBJS)
	$(AR) rcs $(cyln_LIB) $(cyln_OBJS)

