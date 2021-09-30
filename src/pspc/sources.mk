include $(SRC_DIR)/pspc/field/sources.mk
include $(SRC_DIR)/pspc/solvers/sources.mk
include $(SRC_DIR)/pspc/iterator/sources.mk
include $(SRC_DIR)/pspc/sweep/sources.mk

pspc_= \
  $(pspc_field_) \
  $(pspc_solvers_) \
  $(pspc_iterator_) \
  $(pspc_sweep_) \
  pspc/System.cpp 

pspc_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspc_))
pspc_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_:.cpp=.o))

$(pspc_LIB): $(pspc_OBJS)
	$(AR) rcs $(pspc_LIB) $(pspc_OBJS)

