include $(SRC_DIR)/r1d/domain/sources.mk
include $(SRC_DIR)/r1d/solvers/sources.mk
include $(SRC_DIR)/r1d/iterator/sources.mk
include $(SRC_DIR)/r1d/sweep/sources.mk
include $(SRC_DIR)/r1d/misc/sources.mk

r1d_=\
  $(r1d_domain_) \
  $(r1d_solvers_) \
  $(r1d_iterator_) \
  $(r1d_sweep_) \
  $(r1d_misc_) \
  r1d/System.cpp \
  r1d/SystemAccess.cpp 

r1d_OBJS=\
     $(addprefix $(BLD_DIR)/, $(r1d_:.cpp=.o))

$(r1d_LIB): $(r1d_OBJS)
	$(AR) rcs $(r1d_LIB) $(r1d_OBJS)

