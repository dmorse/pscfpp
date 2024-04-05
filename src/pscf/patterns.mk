# ---------------------------------------------------------------------
# File: src/pscf/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at the src/pscf directory, which
# contains files that define classes and functions in the Pscf namespace. 
# It is included by all "makefile" files in this directory tree. 
#-----------------------------------------------------------------------

# List of PSCF-specific libraries needed in src/pscf (the order matters)
PSCF_LIBS=$(pscf_LIB) $(util_LIB) 

# All libraries needed by executables in src/pscf (includes external)
LIBS=$(PSCF_LIBS)

# Add paths to Gnu Scientific Library (GSL)
INCLUDES+=$(GSL_INC)
LIBS+=$(GSL_LIB) 

# Conditionally enable OpenMP
ifdef PSCF_OPENMP
  CXXFLAGS+=$(OPENMP_FLAGS)
  LDFLAGS+=(OPENMP_FLAGS)
  INCLUDES+=$(OPENMP_INC)
  LIBS+=$(OPENMP_LIB) 
endif

# Conditionally add CUDA header include and library paths
ifdef PSCF_CUDA
  INCLUDES+=$(CUDA_INC)
  LIBS+=$(CUDA_LIB)
endif

# Preprocessor macro definitions needed in src/pscf
DEFINES=$(PSCF_DEFS) $(UTIL_DEFS)

# Arguments for MAKEDEP
MAKEDEP_ARGS=$(CPPFLAGS) $(INCLUDES) $(DEFINES)
MAKEDEP_ARGS+= -A$(BLD_DIR)/config.mk
MAKEDEP_ARGS+= -A$(BLD_DIR)/util/config.mk
MAKEDEP_ARGS+= -S$(SRC_DIR)
MAKEDEP_ARGS+= -B$(SRC_DIR)

# Arguments for MAKEDEP for C++
MAKEDEP_CXX_ARGS=$(MAKEDEP_ARGS)
ifdef PSCF_OPENMP
  MAKEDEP_CXX_ARGS+=$(OPENMP_FLAGS)
endif

# Pattern rule to compile *.cpp C++ source files in src/pscf
# Note: Creates a *.d dependency file as a side effect of compilation
$(BLD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
   ifdef MAKEDEP
	$(MAKEDEP) $(MAKEDEP_CMD) $(MAKEDEP_CXX_ARGS) $<
   endif

# Pattern rule to compile *.cu CUDA source files in src/pscf
# Note: Creates a *.d dependency file as a side effect of compilation
$(BLD_DIR)/%.o: $(SRC_DIR)/%.cu
	$(NVXX) $(CPPFLAGS) $(NVXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
   ifdef MAKEDEP_CUDA
	$(MAKEDEP_CUDA) $(MAKEDEP_CUDA_CMD) $(MAKEDEP_ARGS) $<
   endif

# Pattern rule to compile Test programs in src/pscf/tests
$(BLD_DIR)/%Test: $(BLD_DIR)/%Test.o $(PSCF_LIBS)
ifdef PSCF_CUDA
	$(NVXX) $(LDFLAGS) -o $@ $< $(LIBS)
else
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
endif

# Note: In the linking rule for tests, we include the list $(PSCF_LIBS) of 
# PSCF-specific libraries as dependencies but link to the list $(LIBS) 
# that can include external libraries
