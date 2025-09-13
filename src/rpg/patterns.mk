# ---------------------------------------------------------------------
# File: src/rpg/patterns.mk
#
# This makefile fragment contains the pattern rule used to compile all 
# sources files in the directory tree rooted at the src/rpg directory. 
# It is included by all makefiles in this directory tree. 
#-----------------------------------------------------------------------
# Variables used in pattern rules

# List of relevant static libraries defined by PSCF (the order matters)
PSCF_LIBS=$(rpg_LIB) $(prdc_LIB) $(pscf_LIB) $(util_LIB)

# List of all libraries needed in src/rpg (including external libraries)
LIBS=$(PSCF_LIBS)

# Add header include and ibrary paths to Gnu scientific library (GSL)
INCLUDES+=$(GSL_INC)
LIBS+=$(GSL_LIB)

# Add header include and library paths to CUDA FFT library
INCLUDES+=$(CUDA_INC)
LIBS+=$(CUDA_LIB)

# Arguments for MAKEDEP
MAKEDEP_ARGS=$(CPPFLAGS) $(INCLUDES)
MAKEDEP_ARGS+= -A$(BLD_DIR)/config.mk
MAKEDEP_ARGS+= -S$(SRC_DIR)
MAKEDEP_ARGS+= -B$(BLD_DIR)

#-----------------------------------------------------------------------
# Pattern rules

# Pattern rule to compile *.cpp C++ source files in src/rpg
# Note: Creates a *.d dependency file as a side effect 
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	@SDIR=$$(dirname "$@"); if [ ! -d "$$SDIR" ]; then mkdir -p "$$SDIR"; fi
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(CXXFLAGS) -c -o $@ $<
	$(MAKEDEP) $(MAKEDEP_CMD) $(MAKEDEP_ARGS) $<

# Pattern rule to compile *.cu CUDA source files in src/rpg
# Note: Creates a *.d dependency file as a side effect 
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cu
	@SDIR=$$(dirname "$@"); if [ ! -d "$$SDIR" ]; then mkdir -p "$$SDIR"; fi
	$(NVXX) $(CPPFLAGS) $(INCLUDES) $(NVXXFLAGS) -c -o $@ $<
	$(MAKEDEP_CUDA) $(MAKEDEP_CUDA_CMD) $(MAKEDEP_ARGS) $<

# Pattern rule to link executable Test programs in src/rpg/tests
$(BLD_DIR)/%Test: $(BLD_DIR)/%Test.o  $(PSCF_LIBS)
	$(NVXX) $(LDFLAGS) -o $@ $< $(LIBS)

# Note: In the linking rule for unit test programs, we include the list 
# $(PSCF_LIBS) of PSCF-specific libraries as prerequisites but link to 
# the full list $(LIBS) of libraries that includes all relevant external 
# libraries.
