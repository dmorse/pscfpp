# ---------------------------------------------------------------------
# File: src/fd1d/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at the src/fd1d directory, which
# contains all source code for the Util namespace. It is included by
# all "makefile" files in this directory tree. 
#
# This file should be included in other makefiles after inclusion of
# the files src/config.mk and src/fd1d/config.mk because this file
# uses makefile variables defined in those files.
#-----------------------------------------------------------------------

# All libraries needed in src/fd1d
LIBS=$(fd1d_LIB) $(pscf_LIB) $(util_LIB)

# Preprocessor macro definitions needed in src/fd1d
DEFINES=$(PSCF_DEFS) $(UTIL_DEFS)

# Dependencies on build configuration files
MAKE_DEPS= -A$(BLD_DIR)/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/fd1d/config.mk

# Pattern rule to compile *.cpp class source files in src/fd1d
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Pattern rule to compile *.cc test programs in src/fd1d/tests
$(BLD_DIR)/% $(BLD_DIR)/%.o:$(SRC_DIR)/%.cc $(LIBS)
	$(CXX) $(TESTFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $(@:.o=) $@ $(LIBS) $(PSCF_GSL_LIB)
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

