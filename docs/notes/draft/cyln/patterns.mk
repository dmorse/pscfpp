# ---------------------------------------------------------------------
# File: src/cyln/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at the src/cyln directory, which
# contains all source code for the PsSp namespace. It is included by
# all "makefile" files in this directory tree. 
#
# This file must be included in other makefiles after inclusion of
# the root src/config.mk and relevant namespace level config.mk files 
# in the build directory, because this file uses makefile variables 
# defined in those configuration files.
#-----------------------------------------------------------------------

# Local pscf-specific libraries needed in src/cyln
PSCF_LIBS=$(cyln_LIB) $(pscf_LIB) $(util_LIB)

# All libraries needed in executables built in src/cyln
LIBS=$(PSCF_LIBS)
LIBS+=$(GSL_LIB) 
ifdef CYLN_FFTW
LIBS+=$(CYLN_FFTW_LIB) 
endif

# Preprocessor macro definitions needed in src/cyln
DEFINES=$(UTIL_DEFS) $(PSCF_DEFS) $(CYLN_DEFS) 

# Dependencies on build configuration files
MAKE_DEPS= -A$(BLD_DIR)/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/util/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/pscf/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/cyln/config.mk

# Pattern rule to compile *.cpp class source files in src/cyln
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Pattern rule to compile *.cc test programs in src/cyln/tests
$(BLD_DIR)/% $(BLD_DIR)/%.o:$(SRC_DIR)/%.cc $(PSCF_LIBS)
	$(CXX) $(TESTFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $(@:.o=) $@ $(LIBS)
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

