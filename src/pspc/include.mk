# ---------------------------------------------------------------------- #
# This file must be included by every makefile in the pspc/ directory.   #
# It must be included after the config.mk file located in the root of    #
# the build  directory (referred to by a relative path), which defines   #
# values for the macros $(BLD_DIR) and $(SRC_DIR) used below.            # 
# ---------------------------------------------------------------------- #
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/pscf/config.mk
include $(BLD_DIR)/prdc/config.mk
include $(BLD_DIR)/pspc/config.mk
include $(SRC_DIR)/pspc/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/pspc/sources.mk
