# ---------------------------------------------------------------------
# File: src/prdc/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at the src/prdc directory, which
# contains all source code for the Pscf::Prdc namespace. It is included 
# by all "makefile" files in the src/prdc directory tree. 
#-----------------------------------------------------------------------
# Variable definitions used in pattern rules

# Local pscf-specific libraries needed in src/prdc
PRDC_LIBS=$(prdc_LIB) $(pscf_LIB) $(util_LIB)

# All libraries needed by executables in src/prdc (including external)
LIBS=$(PRDC_LIBS)

# Add paths to Gnu scientific library (GSL)
INCLUDES+=$(GSL_INC)
LIBS+=$(GSL_LIB) 

# Add paths to FFTW Fast Fourier transform library
INCLUDES+=$(FFTW_INC)
LIBS+=$(FFTW_LIB) 

# Add paths to CUDA FFT library
ifdef PSCF_CUDA
  INCLUDES+=$(CUDA_INC)
  LIBS+=$(CUDA_LIB)
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
	@SDIR=$$(dirname "$@"); if [!-d"$$SDIR"]; then mkdir -p"$$SDIR"; fi
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<
   ifdef MAKEDEP
	$(MAKEDEP) $(MAKEDEP_CMD) $(MAKEDEP_ARGS) $<
   endif

# Pattern rule to compile *.cu class source files in src/prdc
# Note: Creates a *.d dependency file as a side effect of compilation
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cu
	$(NVXX) $(CPPFLAGS) $(NVXXFLAGS) $(INCLUDES) -c -o $@ $<
   ifdef MAKEDEP_CUDA
	$(MAKEDEP_CUDA) $(MAKEDEP_CUDA_CMD) $(MAKEDEP_ARGS) $<
   endif

# Pattern rule to create executable Test programs in src/prdc/tests
$(BLD_DIR)/%Test: $(BLD_DIR)/%Test.o $(PRDC_LIBS)
ifdef PSCF_CUDA
	$(NVXX) $(LDFLAGS) -o $@ $< $(LIBS)
else
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
endif

# Note: In the linking rule for tests, we include the list $(PRDC_LIBS) 
# of PSCF-specific libraries as dependencies but link to the list $(LIBS) 
# of libraries that includes relevant external libraries
