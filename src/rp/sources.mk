#-----------------------------------------------------------------------
# Source and object file lists for src/rp 

# Include source list files from subdirectories
include $(SRC_DIR)/rp/field/sources.mk
#include $(SRC_DIR)/rp/solvers/sources.mk
#include $(SRC_DIR)/rp/environment/sources.mk
#include $(SRC_DIR)/rp/scft/sources.mk
#include $(SRC_DIR)/rp/fts/sources.mk
#include $(SRC_DIR)/rp/system/sources.mk

# List of compilable C++ source files in src/rp
rp_CPP= \
  $(rp_field_CPP) 
#  $(rp_solvers_CPP) \
#  $(rp_environment_CPP) \
#  $(rp_scft_CPP) \
#  $(rp_fts_CPP) \
#  $(rp_system_CPP) 

# List of compilable CUDA source files in src/rp
rp_CUDA= \
  $(rp_field_CUDA) 
#  $(rp_solvers_CUDA) \
#  $(rp_environment_CUDA) \
#  $(rp_scft_CUDA) \
#  $(rp_fts_CUDA) \
#  $(rp_system_CUDA) 

# List of all object file targets
rp_OBJS=
ifdef PSCF_CPP
   rp_OBJS+= \
       $(addprefix $(BLD_DIR)/, $(rp_CPP:.cpp=.o))
endif
ifdef PSCF_CUDA
   rp_OBJS+= \
       $(addprefix $(BLD_DIR)/, $(rp_CUDA:.cu=.o))
endif

#-----------------------------------------------------------------------
# Path and makefile target for the rp/librp.a library file

rp_LIBNAME=rp
rp_LIB=$(BLD_DIR)/rp/lib$(rp_LIBNAME).a

$(rp_LIB): $(rp_OBJS)
	$(AR) rcs $(rp_LIB) $(rp_OBJS)

#-----------------------------------------------------------------------
