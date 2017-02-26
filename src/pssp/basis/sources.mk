pssp_basis_= \
  pssp/basis/Basis.cpp 

pssp_basis_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_basis_))
pssp_basis_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_basis_:.cpp=.o))

