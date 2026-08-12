# Include source files lists from subdirectories
include $(SRC_DIR)/rp/fts/simulator/sources.mk
include $(SRC_DIR)/rp/fts/compressor/sources.mk
include $(SRC_DIR)/rp/fts/brownian/sources.mk
include $(SRC_DIR)/rp/fts/montecarlo/sources.mk
include $(SRC_DIR)/rp/fts/perturbation/sources.mk
include $(SRC_DIR)/rp/fts/ramp/sources.mk
include $(SRC_DIR)/rp/fts/trajectory/sources.mk
include $(SRC_DIR)/rp/fts/analyzer/sources.mk

# Define object and dependency list variables for rp/fts
rp_fts_OBJS=
rp_fts_DEPS=

rp_fts_OBJS= \
  $(rp_fts_simulator_OBJS) \
  $(rp_fts_compressor_OBJS) \
  $(rp_fts_brownian_OBJS) \
  $(rp_fts_montecarlo_OBJS) \
  $(rp_fts_perturbation_OBJS) \
  $(rp_fts_ramp_OBJS) \
  $(rp_fts_trajectory_OBJS) \
  $(rp_fts_analyzer_OBJS)

rp_fts_DEPS= \
  $(rp_fts_simulator_DEPS) \
  $(rp_fts_compressor_DEPS) \
  $(rp_fts_brownian_DEPS) \
  $(rp_fts_montecarlo_DEPS) \
  $(rp_fts_perturbation_DEPS) \
  $(rp_fts_ramp_DEPS) \
  $(rp_fts_trajectory_DEPS) \
  $(rp_fts_analyzer_DEPS)

