# ----------------------------------------------------------------------- #
# This file must be included by all makefiles in the src/pscf/ directory. #
# It must be included after the config.mk in the root of the build        #
# directory (referred to by an absolute path), which defines values for   #
# the macros $(BLD_DIR) and $(SRC_DIR) as absolute paths.                 #
# ----------------------------------------------------------------------- #
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/pscf/config.mk
include $(SRC_DIR)/pscf/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
