# ---------------------------------------------------------------------- #
# This file must be included by every makefile in the rpc/ directory,    #
# after including the config.mk file in the root of the build directory  #
# ---------------------------------------------------------------------- #
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/rpc/sources.mk
include $(SRC_DIR)/rpc/patterns.mk
