include $(SRC_DIR)/pssp/basis/sources.mk
include $(SRC_DIR)/pssp/field/sources.mk
include $(SRC_DIR)/pssp/iterator/sources.mk
include $(SRC_DIR)/pssp/solvers/sources.mk

pssp_= \
  $(pssp_basis_) \
  $(pssp_field_) \
  $(pssp_solvers_) \
  $(pssp_iterator_) \
  pssp/System.cpp 

pssp_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_))
pssp_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_:.cpp=.o))

$(pssp_LIB): $(pssp_OBJS)
	$(AR) rcs $(pssp_LIB) $(pssp_OBJS)

