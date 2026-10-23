include config.mk

.PHONY: all util pscf r1 rp test \
        clean clean-tests veryclean

# ======================================================================
# Main build targets

# Build all code
all: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd r1; $(MAKE) all
	cd rp; $(MAKE) all

# Build code in Util names (general scientific utilities)
util: 
	cd util; $(MAKE) all

# Build pscf and its dependencies
pscf: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all

# Build pscf_1d program (install in BIN_DIR)
pscf_1d: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd r1; $(MAKE) all

# Build pscf_rpc and/or pscf_rpg (install in BIN_DIR)
rp: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd rp; $(MAKE) all

# ======================================================================
# Test targets

# Run all unit tests
test:
	@cd util/tests; $(MAKE) all; $(MAKE) run-quiet
	@cd pscf/tests; $(MAKE) all; $(MAKE) run-quiet
	@cd r1/tests; $(MAKE) all; $(MAKE) run-quiet
ifdef PSCF_CPP
	@cd rp/tests/cpp; $(MAKE) all; $(MAKE) run-quiet
endif
ifdef PSCF_CUDA
	@cd rp/tests/cuda; $(MAKE) all; $(MAKE) run-quiet
endif
	@cat util/tests/count > count
	@cat pscf/tests/count >> count
	@cat r1/tests/count >> count
ifdef PSCF_CPP
	@cat rp/tests/cpp/count >> count
endif
ifdef PSCF_CUDA
	@cat rp/tests/cuda/count >> count
endif
	@echo " "
	@echo "Summary"
	@cat count
	@rm -f count

# ======================================================================
# Clean targets

# Remove object (*.o), dependency (*.d) and library (*.a) files
clean:
	cd util; $(MAKE) clean
	cd pscf; $(MAKE) clean
	cd r1;  $(MAKE) clean
	cd rp;   $(MAKE) clean

# Clean unit tests
clean-tests:
	cd util/tests; $(MAKE) clean
	cd pscf/tests; $(MAKE) clean
	cd r1/tests; $(MAKE) clean
	cd rp/tests/cpp;  $(MAKE) clean
	cd rp/tests/cuda;  $(MAKE) clean

# Remove all automatically generated files, recreate initial state
veryclean:
	cd util; $(MAKE) veryclean
	cd pscf; $(MAKE) veryclean
	cd r1; $(MAKE) clean
	cd rp;  $(MAKE) clean
	rm -f util/config.mk
ifneq ($(BLD_DIR),$(SRC_DIR))
	rm -f util/makefile
	rm -f pscf/makefile
	rm -f r1/makefile
	rm -f rp/makefile
	rm -f util/tests/makefile
	rm -f pscf/tests/makefile
	rm -f r1/tests/makefile
	rm -f rp/tests/cpp/makefile
	rm -f rp/tests/cuda/makefile
	rm -f setopts
endif
	rm -f config.mk

# ======================================================================
