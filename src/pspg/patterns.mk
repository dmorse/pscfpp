# ---------------------------------------------------------------------
# File: src/pspg/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at the src/pspg directory, which
# contains all source code for the PsSp namespace. It is included by
# all "makefile" files in this directory tree. 
#
# This file must be included in other makefiles after inclusion of
# the root src/config.mk and relevant namespace level config.mk files 
# in the build directory, because this file uses makefile variables 
# defined in those configuration files.
#-----------------------------------------------------------------------

# Local pscf-specific libraries needed in src/pspg (the order matters)
PSPG_LIBS=$(pspg_LIB) $(prdc_LIB) $(pscf_LIB) $(util_LIB)

# All libraries needed in executables built in src/pspg
LIBS=$(PSPG_LIBS)

# Add paths to Gnu scientific library (GSL)
INCLUDES+=$(GSL_INC)
LIBS+=$(GSL_LIB) 

# Add paths to CUDA FFT library
INCLUDES+=$(CUFFT_INC)
LIBS+=$(CUFFT_LIB)

# Preprocessor macro definitions specific to pspg/ directory 
PSPG_DEFS+=-DPSPG_FFTW -DGPU_OUTER

# Preprocessor macro definitions needed in src/pspg
DEFINES=$(UTIL_DEFS) $(PSCF_DEFS) $(PRDC_DEFS) $(PSPG_DEFS) 

# Dependencies on build configuration files
MAKE_DEPS= -A$(BLD_DIR)/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/util/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/pscf/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/prdc/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/pspg/config.mk

# Pattern rule to compile *.cpp C++ source files in src/pspg
# Note: Creates a *.d dependency file as a side effect 
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
   ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
   endif

# Pattern rule to compile *.cu CUDA source files in src/pspg
# Note: Creates a *.d dependency file as a side effect 
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cu
	$(NVXX) $(CPPFLAGS) $(NVXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
   ifdef MAKEDEP_CUDA
	$(MAKEDEP_CUDA) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
   endif

# Pattern rule to link executable Test programs in src/pspg/tests
$(BLD_DIR)/%Test: $(BLD_DIR)/%Test.o  $(PSPG_LIBS)
	$(NVXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $@ $< $(LIBS)
