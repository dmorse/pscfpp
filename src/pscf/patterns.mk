# ---------------------------------------------------------------------
# File: src/pscf/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at the src/pscf directory, which
# contains files that define classes and functions in the Pscf namespace. 
# It is included by all "makefile" files in the src/pscf directory tree. 
#-----------------------------------------------------------------------
# Variable definitions used in pattern rules

# List of PSCF-specific libraries needed in src/pscf (the order matters)
PSCF_LIBS=$(pscf_LIB) $(util_LIB)

# All libraries needed by executables in src/pscf (includes external)
LIBS=$(PSCF_LIBS)

# Add paths to Gnu Scientific Library (GSL)
INCLUDES+=$(GSL_INC)
LIBS+=$(GSL_LIB) 

# Conditionally add CUDA header include and library paths
ifdef PSCF_CUDA
  INCLUDES+=$(CUDA_INC)
  LIBS+=$(CUDA_LIB)
endif

# Arguments for MAKEDEP
MAKEDEP_ARGS=$(CPPFLAGS) $(INCLUDES)
MAKEDEP_ARGS+= -A$(BLD_DIR)/config.mk
MAKEDEP_ARGS+= -S$(SRC_DIR)
MAKEDEP_ARGS+= -B$(BLD_DIR)

# Arguments for MAKEDEP for C++
MAKEDEP_CXX_ARGS=$(MAKEDEP_ARGS)

#-----------------------------------------------------------------------
# Pattern rules

# Pattern rule to compile *.cpp C++ source files in src/pscf
# Note: Creates a *.d dependency file as a side effect of compilation
$(BLD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@SDIR=$$(dirname "$@"); if [ ! -d "$$SDIR" ]; then mkdir -p "$$SDIR"; fi
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<
   ifdef MAKEDEP
	$(MAKEDEP) $(MAKEDEP_CMD) $(MAKEDEP_CXX_ARGS) $<
   endif

# Pattern rule to compile *.cu CUDA source files in src/pscf
# Note: Creates a *.d dependency file as a side effect of compilation
$(BLD_DIR)/%.o: $(SRC_DIR)/%.cu
	@SDIR=$$(dirname "$@"); if [ ! -d "$$SDIR" ]; then mkdir -p "$$SDIR"; fi
	$(NVXX) $(CPPFLAGS) $(NVXXFLAGS) $(INCLUDES) -c -o $@ $<
   ifdef MAKEDEP_CUDA
	$(MAKEDEP_CUDA) $(MAKEDEP_CUDA_CMD) $(MAKEDEP_ARGS) $<
   endif

# Pattern rule to create exectuable Test programs in src/pscf/tests
$(BLD_DIR)/%Test: $(BLD_DIR)/%Test.o $(PSCF_LIBS)
ifdef PSCF_CUDA
	$(NVXX) $(LDFLAGS) -o $@ $< $(LIBS)
else
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
endif

# Note: In the linking rule for tests, we include the list $(PSCF_LIBS) of 
# PSCF-specific libraries as dependencies but link to the list $(LIBS) 
# that can include external libraries
