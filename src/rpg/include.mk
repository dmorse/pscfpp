# ---------------------------------------------------------------------- #
# This file should be included by every makefile in the rpg/ directory.  #
# It must be included after the config.mk file in the root of the build  #
# directory (referred to by a relative path), which defines values for   #
# the macros $(BLD_DIR) and $(SRC_DIR) that are used below.              #
# ---------------------------------------------------------------------- #
include $(BLD_DIR)/util/config.mk

include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/rpg/sources.mk

include $(SRC_DIR)/rpg/patterns.mk
