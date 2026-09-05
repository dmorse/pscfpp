#-----------------------------------------------------------------------
# Source and object file lists for src/pscf 

include $(SRC_DIR)/pscf/backend/sources.mk
include $(SRC_DIR)/pscf/math/sources.mk
include $(SRC_DIR)/pscf/mesh/sources.mk
include $(SRC_DIR)/pscf/chem/sources.mk
include $(SRC_DIR)/pscf/interaction/sources.mk
include $(SRC_DIR)/pscf/floryHuggins/sources.mk
include $(SRC_DIR)/pscf/correlation/sources.mk
include $(SRC_DIR)/pscf/environment/sources.mk
include $(SRC_DIR)/pscf/iterator/sources.mk
include $(SRC_DIR)/pscf/sweep/sources.mk

# Standard C++ files
pscf_OBJS= \
  $(pscf_backend_OBJS) \
  $(pscf_math_OBJS) $(pscf_mesh_OBJS) \
  $(pscf_chem_OBJS) $(pscf_interaction_OBJS) \
  $(pscf_floryHuggins_OBJS) $(pscf_correlation_OBJS) \
  $(pscf_environment_OBJS) $(pscf_iterator_OBJS) \
  $(pscf_sweep_OBJS)

pscf_DEPS= \
  $(pscf_backend_DEPS) \
  $(pscf_math_DEPS) $(pscf_mesh_DEPS) \
  $(pscf_chem_DEPS) $(pscf_interaction_DEPS) \
  $(pscf_floryHuggins_DEPS) $(pscf_correlation_DEPS) \
  $(pscf_environment_DEPS) $(pscf_iterator_DEPS) \
  $(pscf_sweep_DEPS)

#-----------------------------------------------------------------------
# Path and rule for libpscf.a library 

pscf_LIB=$(BLD_DIR)/pscf/libpscf.a

$(pscf_LIB): $(pscf_OBJS)
	$(AR) rcs $(pscf_LIB) $(pscf_OBJS) 

#-----------------------------------------------------------------------
