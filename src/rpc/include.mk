# ---------------------------------------------------------------------- #
# This file must be included by every makefile in the rpc/ directory.   #
# It must be included after the config.mk file located in the root of    #
# the build  directory (referred to by a relative path), which defines   #
# values for the macros $(BLD_DIR) and $(SRC_DIR) used below.            # 
# ---------------------------------------------------------------------- #
include $(BLD_DIR)/util/config.mk

include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/rpc/sources.mk

include $(SRC_DIR)/rpc/patterns.mk
