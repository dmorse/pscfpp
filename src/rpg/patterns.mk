# ---------------------------------------------------------------------
# File: src/rpg/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at the src/rpg directory, which
# contains all source code for the PsSp namespace. It is included by
# all "makefile" files in this directory tree. 
#
# This file must be included in other makefiles after inclusion of
# the root src/config.mk and relevant namespace level config.mk files 
# in the build directory, because this file uses makefile variables 
# defined in those configuration files.
#-----------------------------------------------------------------------

# Local PSCF-specific libraries needed in src/rpg (the order matters)
PSCF_LIBS=$(rpg_LIB) $(prdc_LIB) $(pscf_LIB) $(util_LIB)

# All libraries needed by executables in src/rpg (including external libs)
LIBS=$(PSCF_LIBS)

# Add paths to Gnu scientific library (GSL)
INCLUDES+=$(GSL_INC)
LIBS+=$(GSL_LIB) 

# Add paths to CUDA FFT library
INCLUDES+=$(CUDA_INC)
LIBS+=$(CUDA_LIB)

# Preprocessor macro definitions specific to rpg/ directory 
RPG_DEFS+=-DRPG_FFTW -DGPU_OUTER

# Preprocessor macro definitions needed in src/rpg
DEFINES=$(UTIL_DEFS) $(PSCF_DEFS) $(PRDC_DEFS) $(RPG_DEFS) 

# Arguments for MAKEDEP
MAKEDEP_ARGS=$(CPPFLAGS) $(INCLUDES) $(DEFINES)
MAKEDEP_ARGS+= -A$(BLD_DIR)/config.mk
MAKEDEP_ARGS+= -A$(BLD_DIR)/util/config.mk
MAKEDEP_ARGS+= -A$(BLD_DIR)/pscf/config.mk
MAKEDEP_ARGS+= -A$(BLD_DIR)/prdc/config.mk
MAKEDEP_ARGS+= -A$(BLD_DIR)/rpg/config.mk
MAKEDEP_ARGS+= -S$(SRC_DIR)
MAKEDEP_ARGS+= -B$(SRC_DIR)

# Pattern rule to compile *.cpp C++ source files in src/rpg
# Note: Creates a *.d dependency file as a side effect 
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
   ifdef MAKEDEP
	$(MAKEDEP) $(MAKEDEP_CMD) $(MAKEDEP_ARGS) $<
   endif

# Pattern rule to compile *.cu CUDA source files in src/rpg
# Note: Creates a *.d dependency file as a side effect 
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cu
	$(NVXX) $(CPPFLAGS) $(NVXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
   ifdef MAKEDEP_CUDA
	$(MAKEDEP_CUDA) $(MAKEDEP_CUDA_CMD) $(MAKEDEP_ARGS) $<
   endif

# Pattern rule to link executable Test programs in src/rpg/tests
$(BLD_DIR)/%Test: $(BLD_DIR)/%Test.o  $(PSCF_LIBS)
	$(NVXX) $(LDFLAGS) -o $@ $< $(LIBS)

# Note: In the linking rule for tests, we include the list $(PSCF_LIBS) 
# of PSCF-specific libraries as dependencies but link to the list 
# $(LIBS) of libraries that includes any relevant external libraries.