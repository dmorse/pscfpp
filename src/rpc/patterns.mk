# ---------------------------------------------------------------------
# File: src/rpc/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at the src/rpc directory, which
# contains source files defined in the Pscf::Rpc namespace. It is 
# included by all makefiles in this directory tree. 
#-----------------------------------------------------------------------

# PSCF-specific static libraries needed in src/rpc (the order matters)
# Variables $(rpc_LIB) etc. are defined in namespace config.mk files
PSCF_LIBS= $(rpc_LIB) $(prdc_LIB) $(pscf_LIB) $(util_LIB) 

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
  LDFLAGS+=$(OPENMP_FLAGS)
  INCLUDES+=$(OPENMP_INC)
  LIBS+=$(OPENMP_LIB) 
endif

# List of all preprocessor macro definitions needed in src/rpc
DEFINES=$(UTIL_DEFS) $(PSCF_DEFS) 

# Arguments for MAKEDEP
MAKEDEP_ARGS=$(CPPFLAGS) $(INCLUDES) $(DEFINES)
MAKEDEP_ARGS+= -A$(BLD_DIR)/config.mk
MAKEDEP_ARGS+= -A$(BLD_DIR)/util/config.mk
MAKEDEP_ARGS+= -S$(SRC_DIR)
MAKEDEP_ARGS+= -B$(BLD_DIR)

# Arguments for MAKEDEP for C++
MAKEDEP_CXX_ARGS=$(MAKEDEP_ARGS)
ifdef PSCF_OPENMP
  MAKEDEP_CXX_ARGS+=$(OPENMP_FLAGS)
endif

# Pattern rule to compile *.cpp class source files in src/rpc
# Note: Creates a *.d dependency file as a side effect
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
   ifdef MAKEDEP
	$(MAKEDEP) $(MAKEDEP_CMD) $(MAKEDEP_CXX_ARGS) $<
   endif

# Pattern rule to link executable Test programs in src/rpc/tests
$(BLD_DIR)/%Test: $(BLD_DIR)/%Test.o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)

# Note: In the linking rule for tests, we include the list $(PSCF_LIBS) 
# of PSCF-specific libraries as dependencies but link to the list $(LIBS) 
# that can include external libraries
