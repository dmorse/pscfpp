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

# Local pscf-specific libraries needed in src/pspg
PSPG_LIBS=$(pspg_LIB) $(pscf_LIB) $(util_LIB)

# All libraries needed in executables built in src/pssp
LIBS=$(PSPG_LIBS)

# Add paths to Gnu scientific library (GSL)
INCLUDES+=$(GSL_INC)
LIBS+=$(GSL_LIB) 

# Add paths to Cuda FFT library
PSPG_DEFS+=-DPSPG_FFTW -DGPU_OUTER
INCLUDES+=$(CUFFT_INC)
LIBS+=$(CUFFT_LIB)

# Preprocessor macro definitions needed in src/pssp
DEFINES=$(UTIL_DEFS) $(PSCF_DEFS) $(PSPG_DEFS) 

# Dependencies on build configuration files
MAKE_DEPS= -A$(BLD_DIR)/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/util/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/pscf/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/pspg/config.mk

# Pattern rule to compile *.cpp class source files in src/pssp
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Pattern rule to compile *.cu class source files in src/pssp
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cu
	$(NVXX) $(CPPFLAGS) $(NVXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
#ifdef MAKEDEP
#	$(MAKEDEP_CUDA) $(NVXXFLAGS) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
#endif

# Pattern rule to compile *.ccu test programs in src/pssp/tests
$(BLD_DIR)/% $(BLD_DIR)/%.o:$(SRC_DIR)/%.ccu $(PSPG_LIBS)
	$(NVXX)  $(NVXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
	$(NVXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $(@:.o=) $@ $(LIBS)
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

