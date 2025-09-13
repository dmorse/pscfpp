# ---------------------------------------------------------------------
# File: src/r1d/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at the src/r1d directory, which
# contains all source code for the Pscf::R1d namespace. It is included
# by all "makefile" files in this directory tree. 
#-----------------------------------------------------------------------
# Definitions of variables used in pattern rules

# PSCF-specific libraries needed in src/r1d
PSCF_LIBS=$(r1d_LIB) $(pscf_LIB) $(util_LIB)

# All libraries needed by executables in src/r1d (including external)
LIBS=$(PSCF_LIBS)

# Add paths to Gnu scientific library (GSL)
INCLUDES+=$(GSL_INC)
LIBS+=$(GSL_LIB) 

# Arguments for MAKEDEP script
MAKEDEP_ARGS=$(CPPFLAGS) $(INCLUDES)
MAKEDEP_ARGS+= -A$(BLD_DIR)/config.mk
MAKEDEP_ARGS+= -S$(SRC_DIR)
MAKEDEP_ARGS+= -B$(BLD_DIR)

#-----------------------------------------------------------------------
# Pattern rules

# Pattern rule to compile *.cpp class source files in src/r1d
# Note: Creates a *.d dependency file as a side effect
$(BLD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@SDIR=$$(dirname "$@"); if [ ! -d"$$SDIR" ]; then mkdir -p"$$SDIR"; fi
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(CXXFLAGS) -c -o $@ $<
	$(MAKEDEP) $(MAKEDEP_CMD) $(MAKEDEP_ARGS) $<

# Pattern rule to link Test programs in src/r1d/tests
$(BLD_DIR)/%Test: $(BLD_DIR)/%Test.o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)

# Note: In the linking rule for tests, we include the list $(PSCF_LIBS) 
# of PSCF-specific libraries as dependencies but link to list $(LIBS) of
# libraries that includes relevant external libraries
