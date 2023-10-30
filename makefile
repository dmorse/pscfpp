include src/config.mk
# =========================================================================
.PHONY:  all-cpu  all fd1d pspc pspg \
         html html-dev \
         clean clean-tests clean-bin clean-html veryclean 

# =========================================================================
# Main targets 

# All main targets in this file perform out-of-source compilation in which
# intermediate files are created in the bld/ directory. Each of these 
# targets simply cd's to bld/ directory and then invokes make with the 
# same target name from there. In-source compilation may be performed by 
# instead invoking make with the same target names from within the src/ 
# directory.

# Compile all code, including CUDA code if possible
all:
	cd bld; $(MAKE) all

# Build all programs that runs on a conventional cpu 
all-cpu:
	cd bld; $(MAKE) all-cpu

# Build the pscf_1d 1D finite element SCFT program
pscf_1d:
	cd bld; $(MAKE) pscf_1d

# Build the pscf_pc cpu program for periodic structures 
pscf_pc:
	cd bld; $(MAKE) pscf_pc

ifdef PSCF_CUDA
# Build the pscf_pg gpu program for periodic structures
pscf_pg:
	cd bld; $(MAKE) pscf_pg
endif


# =========================================================================
# HTML Documentation
 
html:
	cd docs; $(MAKE) html

html-dev:
	cd docs; $(MAKE) html-dev

# =========================================================================
# Clean-up targets
# ----------------
# Clean and veryclean targets defined here clean up both the bld/ and src/ 
# directories. The veryclean target also removes most other files that were 
# not in the repository, including makefile created by the setup script.


# Remove *.o, *.d and *.a files in both bld/ and src/ directories
clean:
	cd bld; $(MAKE) clean
	cd src; $(MAKE) clean

# Remove files created by compiling and running unit tests.
clean-tests:
	cd src/; $(MAKE) clean-tests
	cd bld/; $(MAKE) clean-tests

# Remove all executable files in $(BIN_DIR)
clean-bin:
	rm -f $(BIN_DIR)/pscf*

# Remove user-generated documentation files
clean-docs:
	cd docs; $(MAKE) clean

# Remove almost all files that were not in the repository
veryclean:
	make clean-bin
	rm -f $(BIN_DIR)/makeDep*
	make clean-docs
	cd bld; $(MAKE) veryclean
	rm bld/makefile
	cd src; $(MAKE) veryclean
	cd examples; ./clean
	cd lib/python/pscfpp; rm -f *.pyc; rm -Rf __pycache__

# ==========================================================================
