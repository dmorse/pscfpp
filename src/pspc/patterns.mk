# ---------------------------------------------------------------------
# File: src/pspc/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at the src/pspc directory, which
# contains source files defined in the Pscf namespace. It is included 
# by all makefiles in this directory tree. 
#
# This file must be included in other makefiles after inclusion of
# the root src/config.mk and relevant namespace level config.mk files 
# in the build directory, because this file uses makefile variables 
# defined in those configuration files.
#-----------------------------------------------------------------------

# PSCF-specific static libraries needed in src/pspc (the order matters)
# Variables $(pspc_LIB) etc. are defined in namespace config.mk files
PSCF_LIBS= $(pspc_LIB) $(prdc_LIB) $(pscf_LIB) $(util_LIB) 

# All libraries needed by main program pscf_pc (including external libs)
LIBS=$(PSCF_LIBS)

# Add paths to Gnu scientific library
INCLUDES+=$(GSL_INC)
LIBS+=$(GSL_LIB) 

# Add paths to FFTW Fast Fourier transform library
INCLUDES+=$(FFTW_INC)
LIBS+=$(FFTW_LIB) 

# Conditionally enable OpenMP
ifdef PSCF_OPENMP
  CXXFLAGS+=$(OPENMP_FLAGS)
  INCLUDES+=$(OPENMP_INC)
  LIBS+=$(OPENMP_LIB) 
endif

# List of all preprocessor macro definitions needed in src/pspc
# Variables $(PSPC_DEFS) etc are initialized in namespace config.mk files
DEFINES=$(UTIL_DEFS) $(PSCF_DEFS) $(PSPC_DEFS) 

# Dependencies on build configuration files
MAKE_DEPS= -A$(BLD_DIR)/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/util/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/pscf/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/prdc/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/pspc/config.mk

# Pattern rule to compile *.cpp class source files in src/pspc
# Note: Creates a *.d dependency file as a side effect
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
   ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
   endif

# Pattern rule to link executable Test programs in src/pscf/tests
# Note: Only list PSCF-specific libraries as dependencies, but link all
$(BLD_DIR)/%Test: $(BLD_DIR)/%Test.o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $@ $< $(LIBS)

