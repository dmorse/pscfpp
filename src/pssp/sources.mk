include $(SRC_DIR)/pssp/basis/sources.mk
include $(SRC_DIR)/pssp/field/sources.mk

pssp_= \
  $(pssp_basis_) \
  $(pssp_field_) 

pssp_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_))
pssp_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_:.cpp=.o))

$(pssp_LIB): $(pssp_OBJS)
	$(AR) rcs $(pssp_LIB) $(pssp_OBJS)

