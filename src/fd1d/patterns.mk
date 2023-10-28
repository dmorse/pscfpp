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

# PSCF-specific libraries needed in src/fd1d
PSCF_LIBS=$(fd1d_LIB) $(pscf_LIB) $(util_LIB)

# All libraries needed by executables in src/fd1d (including external)
LIBS=$(PSCF_LIBS)

# Add paths to Gnu scientific library (GSL)
INCLUDES+=$(GSL_INC)
LIBS+=$(GSL_LIB) 

# Preprocessor macro definitions needed in src/fd1d
DEFINES=$(PSCF_DEFS) $(UTIL_DEFS)

# Arguments for MAKEDEP script
MAKEDEP_ARGS=$(CPPFLAGS) $(INCLUDES) $(DEFINES)
MAKEDEP_ARGS+= -A$(BLD_DIR)/config.mk
MAKEDEP_ARGS+= -A$(BLD_DIR)/util/config.mk
MAKEDEP_ARGS+= -A$(BLD_DIR)/pscf/config.mk
MAKEDEP_ARGS+= -A$(BLD_DIR)/fd1d/config.mk
MAKEDEP_ARGS+= -S$(SRC_DIR)
MAKEDEP_ARGS+= -B$(SRC_DIR)

# Pattern rule to compile *.cpp class source files in src/fd1d
# Note: Creates a *.d dependency file as a side effect
$(BLD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
   ifdef MAKEDEP
	$(MAKEDEP) $(MAKEDEP_CMD) $(MAKEDEP_ARGS) $<
   endif

# Pattern rule to link Test programs in src/fd1d/tests
$(BLD_DIR)/%Test: $(BLD_DIR)/%Test.o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)

# Note: In the linking rule for tests, we include the list $(PSCF_LIBS) 
# of PSCF-specific libraries as dependencies but link to list $(LIBS) of
# libraries that includes relevant external libraries
