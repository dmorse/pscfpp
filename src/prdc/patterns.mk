# ---------------------------------------------------------------------
# File: src/prdc/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at the src/prdc directory, which
# contains all source code for the Pscf::Prdc namespace. It is included 
# by all "makefile" files in the src/prdc directory tree. 
#-----------------------------------------------------------------------
# Variable definitions used in pattern rules

# Local PSCF static libraries needed in src/prdc
# Variables *_LIB give library paths, are defined in sources.mk files
PRDC_LIBS=$(prdc_LIB) $(pscf_LIB) $(util_LIB)

# All libraries needed by executables in src/prdc (including external)
LIBS=$(PRDC_LIBS)

# Add header include and library paths to Gnu scientific library (GSL)
INCLUDES+=$(GSL_INC)
LIBS+=$(GSL_LIB) 

# Add header include, library paths to FFTW Fast Fourier transform library
INCLUDES+=$(FFTW_INC)
LIBS+=$(FFTW_LIB) 

ifdef PSCF_CUDA
  # Add include and library paths to CUDA FFT library
  INCLUDES+=$(CUDA_INC)
  LIBS+=$(CUDA_LIB)
endif

# Define command CXXTEST used to create test executables
ifdef PSCF_CUDA
  CXXTEST = $(NVXX)
else
  CXXTEST = $(CXX)
endif

# Arguments for MAKEDEP
MAKEDEP_ARGS=$(CPPFLAGS) $(INCLUDES) 
MAKEDEP_ARGS+= -A$(BLD_DIR)/config.mk
MAKEDEP_ARGS+= -S$(SRC_DIR)
MAKEDEP_ARGS+= -B$(BLD_DIR)

#-----------------------------------------------------------------------
# Pattern rules

# Pattern rule to compile *.cpp class source files in src/prdc
# Note: Creates a *.d dependency file as a side effect of compilation
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	@SDIR=$$(dirname "$@"); if [ ! -d "$$SDIR" ]; then mkdir -p "$$SDIR"; fi
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(CXXFLAGS) -c -o $@ $<
	$(MAKEDEP) $(MAKEDEP_CMD) $(MAKEDEP_ARGS) $<

# Pattern rule to compile *.cu class source files in src/prdc
# Note: Creates a *.d dependency file as a side effect of compilation
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cu
	@SDIR=$$(dirname "$@"); if [ ! -d "$$SDIR" ]; then mkdir -p "$$SDIR"; fi
	$(NVXX) $(CPPFLAGS) $(INCLUDES) $(NVXXFLAGS) -c -o $@ $<
	$(MAKEDEP_CUDA) $(MAKEDEP_CUDA_CMD) $(MAKEDEP_ARGS) $<

# Pattern rule to create executable Test programs in src/prdc/tests
$(BLD_DIR)/%Test: $(BLD_DIR)/%Test.o $(PRDC_LIBS)
	$(CXXTEST) $(LDFLAGS) -o $@ $< $(LIBS)

# Note: We use $(LIBS) in the recipe for unit test executables but
# $(PSCF_LIBS) in the prerequisite list, so that all libraries are
# linked but the build system is only responsible for updating those
# that are part of PSCF.
