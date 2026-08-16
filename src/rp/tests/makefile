include ../../config.mk

.PHONY: all run run-quite clean-outputs clean

# ======================================================================
# Main build targets

# Compile unit tests
all: 
ifdef PSCF_CPP
	cd cpp; $(MAKE) all
endif
ifdef PSCF_CUDA
	cd cuda; $(MAKE) all
endif

run: 
ifdef PSCF_CPP
	cd cpp; $(MAKE) run
endif
ifdef PSCF_CUDA
	cd cuda; $(MAKE) run
endif

run-quiet: 
ifdef PSCF_CPP
	cd cpp; $(MAKE) run-quiet
endif
ifdef PSCF_CUDA
	cd cuda; $(MAKE) run-quiet
endif

clean-outputs: 
	cd cpp; $(MAKE) clean-outputs
	cd cuda; $(MAKE) clean-outputs

clean: 
	cd cpp; $(MAKE) clean
	cd cuda; $(MAKE) clean

