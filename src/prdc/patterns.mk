# ---------------------------------------------------------------------
# File: src/prdc/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at the src/prdc directory, which
# contains all source code for the Pscf::Prdc namespace. It is included 
# by all "makefile" files in this directory tree. 
#
# This file must be included in other makefiles after inclusion of
# the root src/config.mk and relevant namespace level config.mk files 
# in the build directory, because this file uses makefile variables 
# defined in those configuration files.
#-----------------------------------------------------------------------

# Local pscf-specific libraries needed in src/prdc
PRDC_LIBS=$(prdc_LIB) $(pscf_LIB) $(util_LIB)

# All libraries needed in executables built in src/prdc
LIBS=$(PRDC_LIBS)

# Add paths to Gnu scientific library (GSL)
INCLUDES+=$(GSL_INC)
LIBS+=$(GSL_LIB) 

# Add paths to FFTW Fast Fourier transform library
INCLUDES+=$(FFTW_INC)
LIBS+=$(FFTW_LIB) 

# Add paths to CUDA FFT library
ifdef PSCF_CUDA
  PRDC_DEFS+=-DPRDC_FFTW -DGPU_OUTER
  PRDC_CUFFT_LIB=-lcufft -lcudart -lcuda -lcurand
  INCLUDES+=$(CUFFT_INC)
  LIBS+=$(PRDC_CUFFT_LIB)
endif

# Preprocessor macro definitions needed in src/prdc
DEFINES=$(UTIL_DEFS) $(PSCF_DEFS) $(PRDC_DEFS) 

# Dependencies on build configuration files
MAKE_DEPS= -A$(BLD_DIR)/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/util/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/pscf/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/prdc/config.mk

# Pattern rule to compile *.cpp class source files in src/prdc
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Pattern rule to compile *.cu class source files in src/prdc
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cu
ifdef PSCF_CUDA
	$(NVXX) $(CPPFLAGS) $(NVXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
  ifdef MAKEDEP_CUDA
	$(MAKEDEP_CUDA) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
  endif
else
  # Attempt to compile *.cu files as *.cpp files ifndef PSCF_CUDA
  # Make temporary copy of *.cu file with *.cpp extension, compile as C++ file
	cp $< $(TMP) $(<:.cu=.cpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $(<:.cu=.cpp)
  ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $(<:.cu=.cpp)
  endif
  # Remove temporary *.cpp file
	rm $(<:.cu=.cpp)
endif

# Pattern rule to compile Test programs in src/prdc/tests
$(BLD_DIR)/%Test: $(BLD_DIR)/%Test.o $(PRDC_LIBS)
ifdef PSCF_CUDA
	$(NVXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $@ $< $(LIBS)
else
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $@ $< $(LIBS)
endif
