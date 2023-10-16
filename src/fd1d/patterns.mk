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

# Local pscf-specific libraries needed in src/fd1d
FD1D_LIBS=$(fd1d_LIB) $(pscf_LIB) $(util_LIB)

# List of all libraries needed by executables in src/fd1d
LIBS=$(FD1D_LIBS)

# Add paths to Gnu scientific library (GSL)
INCLUDES+=$(GSL_INC)
LIBS+=$(GSL_LIB) 

# Preprocessor macro definitions needed in src/fd1d
DEFINES=$(PSCF_DEFS) $(UTIL_DEFS)

# Dependencies on build configuration files
MAKE_DEPS= -A$(BLD_DIR)/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/util/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/pscf/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/fd1d/config.mk

# Pattern rule to compile *.cpp class source files in src/fd1d
# Note: Creates a *.d dependency file as a side effect
$(BLD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
   ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
   endif

# Pattern rule to link Test programs in src/pscf/tests
$(BLD_DIR)/%Test: $(BLD_DIR)/%Test.o $(FD1D_LIBS)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $@ $< $(LIBS)

