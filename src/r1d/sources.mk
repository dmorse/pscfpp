#-----------------------------------------------------------------------
# Source and object file lists for src/r1d 

# Include source list files from subdirectories
include $(SRC_DIR)/r1d/domain/sources.mk
include $(SRC_DIR)/r1d/solvers/sources.mk
include $(SRC_DIR)/r1d/iterator/sources.mk
include $(SRC_DIR)/r1d/sweep/sources.mk
include $(SRC_DIR)/r1d/misc/sources.mk

# List of source files in src/r1d
r1d_=\
  $(r1d_domain_) \
  $(r1d_solvers_) \
  $(r1d_iterator_) \
  $(r1d_sweep_) \
  $(r1d_misc_) \
  r1d/System.cpp \
  r1d/SystemAccess.cpp 

# List of object file targets
r1d_OBJS=\
     $(addprefix $(BLD_DIR)/, $(r1d_:.cpp=.o))

#-----------------------------------------------------------------------
# Path and target rule for the r1d/libr1d.a library 

r1d_LIBNAME=r1d
r1d_LIB=$(BLD_DIR)/r1d/lib$(r1d_LIBNAME).a

$(r1d_LIB): $(r1d_OBJS)
	$(AR) rcs $(r1d_LIB) $(r1d_OBJS)

#-----------------------------------------------------------------------
