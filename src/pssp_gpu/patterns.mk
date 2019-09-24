# ---------------------------------------------------------------------
# File: src/pssp/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at the src/pssp directory, which
# contains all source code for the PsSp namespace. It is included by
# all "makefile" files in this directory tree. 
#
# This file must be included in other makefiles after inclusion of
# the root src/config.mk and relevant namespace level config.mk files 
# in the build directory, because this file uses makefile variables 
# defined in those configuration files.
#-----------------------------------------------------------------------

# Local pscf-specific libraries needed in src/pssp_gpu
PSCF_LIBS=$(pssp_gpu_LIB) $(pscf_LIB) $(util_LIB)

# All libraries needed in executables built in src/pssp
LIBS=$(PSCF_LIBS)
ifdef PSCF_GSL
LIBS+=$(PSCF_GSL_LIB) 
endif
ifdef PSSP_CUDA
LIBS+=$(PSSP_CUFFT_LIB)
endif

# Preprocessor macro definitions needed in src/pssp
DEFINES=$(UTIL_DEFS) $(PSCF_DEFS) $(PSSP_GPU_DEFS) 

# Dependencies on build configuration files
MAKE_DEPS= -A$(BLD_DIR)/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/util/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/pscf/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/pssp_gpu/config.mk

# Pattern rule to compile *.cpp class source files in src/pssp
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Pattern rule to compile *.cu class source files in src/pssp
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cu
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
#ifdef MAKEDEP
#	$(MAKEDEP) -arch=sm_30 $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
#endif

# Pattern rule to compile *.cc test programs in src/pssp/tests
$(BLD_DIR)/% $(BLD_DIR)/%.o:$(SRC_DIR)/%.cu $(PSCF_LIBS)
	nvcc  -O3 $(INCLUDES) $(DEFINES) -c -o $@ $<
	nvcc $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $(@:.o=) $@ $(LIBS)
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

