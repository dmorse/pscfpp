# ---------------------------------------------------------------------
# File: src/rpc/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at the src/rpc directory. It is
# included by all makefiles in this directory tree. 
#-----------------------------------------------------------------------
# Variable definitions used in pattern rules

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

# Arguments for MAKEDEP
MAKEDEP_ARGS=$(CPPFLAGS) $(INCLUDES)
MAKEDEP_ARGS+= -A$(BLD_DIR)/config.mk
MAKEDEP_ARGS+= -S$(SRC_DIR)
MAKEDEP_ARGS+= -B$(BLD_DIR)

#-----------------------------------------------------------------------
# Pattern rules

# Pattern rule to compile *.cpp class source files in src/rpc
# Note: Creates a *.d dependency file as a side effect
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	@SDIR=$$(dirname "$@"); if [ ! -d "$$SDIR" ]; then mkdir -p "$$SDIR"; fi
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<
   ifdef MAKEDEP
	$(MAKEDEP) $(MAKEDEP_CMD) $(MAKEDEP_ARGS) $<
   endif

# Pattern rule to link executable Test programs in src/rpc/tests
$(BLD_DIR)/%Test: $(BLD_DIR)/%Test.o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)

# Note: In the linking rule for tests, we include the list $(PSCF_LIBS) 
# of PSCF-specific libraries as prerequisites but link to the list 
# $(LIBS) of libraries that includes external libraries

# Note: There are no *.cu CUDA source files in src/rpc, and so no rule to 
# compile them is defined in this directory.

