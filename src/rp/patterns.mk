# ---------------------------------------------------------------------
# File: src/rp/patterns.mk
#
# This makefile fragment contains the pattern rule used to compile all 
# sources files in the directory tree rooted at the src/rp directory. 
# It is included by all makefiles in this directory tree. 
#-----------------------------------------------------------------------
# Variables used in pattern rules

# List of relevant static libraries defined by PSCF (the order matters)
PSCF_LIBS=$(rp_LIB) $(prdc_LIB) $(pscf_LIB) $(util_LIB)

# List of all libraries needed in src/rp (including external libraries)
LIBS=$(PSCF_LIBS)

# Add header include and ibrary paths to Gnu scientific library (GSL)
INCLUDES+=$(GSL_INC)
LIBS+=$(GSL_LIB)

# Conditionally add paths for the C++ FFTW library
ifdef PSCF_CPP
INCLUDES+=$(FFTW_INC)
LIBS+=$(FFTW_LIB)
endif

# Conditionally add paths for the CUDA FFT library
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

# Pattern rule to compile *.cpp C++ source files in src/rp
# Note: Creates a *.d dependency file as a side effect 
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	@SDIR=$$(dirname "$@"); if [ ! -d "$$SDIR" ]; then mkdir -p "$$SDIR"; fi
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(CXXFLAGS) -c -o $@ $<
	$(MAKEDEP) $(MAKEDEP_CMD) $(MAKEDEP_ARGS) $<

# Pattern rule to compile *.cu CUDA source files in src/rp
# Note: Creates a *.d dependency file as a side effect 
$(BLD_DIR)/%.ou:$(SRC_DIR)/%.cu
	@SDIR=$$(dirname "$@"); if [ ! -d "$$SDIR" ]; then mkdir -p "$$SDIR"; fi
	$(NVXX) $(CPPFLAGS) $(INCLUDES) $(NVXXFLAGS) -c -o $@ $<
	$(MAKEDEP_CUDA) $(MAKEDEP_CUDA_CMD) $(MAKEDEP_ARGS) $<

# Pattern rule to link executable Test programs in src/rp/tests
$(BLD_DIR)/%Test: $(BLD_DIR)/%Test.o  $(PSCF_LIBS)
	$(NVXX) $(LDFLAGS) -o $@ $< $(LIBS)

# Note: In the linking rule for unit test programs, the prerequisite list
# contains the list $(PSCF_LIBS) of PSCF-specific libraries, but the action
# links to the full list $(LIBS) of libraries that includes all relevant 
# external libraries.
