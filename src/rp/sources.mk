#-----------------------------------------------------------------------
# Source and object file lists for src/rp 

# Include source list files from subdirectories
include $(SRC_DIR)/rp/field/sources.mk
include $(SRC_DIR)/rp/solvers/sources.mk
include $(SRC_DIR)/rp/system/sources.mk
#include $(SRC_DIR)/rp/environment/sources.mk
#include $(SRC_DIR)/rp/scft/sources.mk
#include $(SRC_DIR)/rp/fts/sources.mk

rp_OBJS= \
     $(rp_field_OBJS) \
     $(rp_solvers_OBJS) \
     $(rp_system_OBJS) 
   #  $(rp_environment_OBJS) \
   #  $(rp_scft_OBJS) \
   #  $(rp_fts_OBJS) \


#-----------------------------------------------------------------------
# Path and makefile target for the rp/librp.a library file

rp_LIBNAME=rp
rp_LIB=$(BLD_DIR)/rp/lib$(rp_LIBNAME).a

$(rp_LIB): $(rp_OBJS)
	$(AR) rcs $(rp_LIB) $(rp_OBJS)

#-----------------------------------------------------------------------
