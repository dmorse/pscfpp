include config.mk

.PHONY: all all-old util pscf prdc r1d rp rpc rpg cpc \
        test test-cpu test-gpu \
        clean clean-tests veryclean
# ======================================================================
# Main build targets

# Build all code
all: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd r1d; $(MAKE) all
	cd prdc; $(MAKE) all
	cd rp; $(MAKE) all

# Build code in Util names (general scientific utilities)
util: 
	cd util; $(MAKE) all

# Build code required by pscf
pscf: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all

# Build code required by prdc
prdc: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd prdc; $(MAKE) all

# Build pscf_r1d 1D finite difference SCFT program (install in BIN_DIR)
pscf_r1d: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd r1d; $(MAKE) all

# Build pscf_rpc and/or pscf_rpg in rp (install in BIN_DIR)
rp: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd prdc; $(MAKE) all
	cd rp; $(MAKE) all

# Build pscf_rpc CPU code for real periodic fields (install in BIN_DIR)
pscf_rpc: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd prdc; $(MAKE) all
	cd rpc; $(MAKE) all

ifdef PSCF_CUDA
# Build pscf_rpc GPU code for real periodic fields (install in BIN_DIR)
pscf_rpg:
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd prdc; $(MAKE) all
	cd rpg; $(MAKE) all
endif

# Build pscf_cpc CPU code for real periodic fields (install in BIN_DIR)
pscf_cpc: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd prdc; $(MAKE) all
	cd cpc; $(MAKE) all

# ======================================================================
# Test targets

# Run all unit tests
test:
	@cd util/tests; $(MAKE) all; $(MAKE) run-quiet
	@cd pscf/tests; $(MAKE) all; $(MAKE) run-quiet
	@cd prdc/tests; $(MAKE) all; $(MAKE) run-quiet
	@cd r1d/tests; $(MAKE) all; $(MAKE) run-quiet
ifdef PSCF_CPP
	@cd rp/tests/cpp; $(MAKE) all; $(MAKE) run-quiet
endif
ifdef PSCF_CUDA
	@cd rp/tests/cuda; $(MAKE) all; $(MAKE) run-quiet
endif
	@cat util/tests/count > count
	@cat pscf/tests/count >> count
	@cat prdc/tests/count >> count
	@cat r1d/tests/count >> count
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

# Run all CPU unit tests
test-cpu:
	@cd util/tests; $(MAKE) all; $(MAKE) run-quiet
	@cd pscf/tests; $(MAKE) all; $(MAKE) run-quiet
	@cd prdc/tests; $(MAKE) all; $(MAKE) run-quiet
	@cd r1d/tests; $(MAKE) all; $(MAKE) run-quiet
ifdef PSCF_CPP
	@cd rp/tests/cpp; $(MAKE) all; $(MAKE) run-quiet
endif
	@cat util/tests/count > count
	@cat pscf/tests/count >> count
	@cat prdc/tests/count >> count
	@cat r1d/tests/count >> count
ifdef PSCF_CPP
	@cat rp/tests/cpp/count >> count
endif
	@echo " "
	@echo "Summary"
	@cat count
	@rm -f count

test-gpu:
	@cd rp/tests/cuda; @$(MAKE) all; $(MAKE) run-quiet
	#@cd rpg/tests; @$(MAKE) all; $(MAKE) run-quiet

# ======================================================================
# Clean targets

# Remove object (*.o), dependency (*.d) and library (*.a) files
clean:
	cd util; $(MAKE) clean
	cd pscf; $(MAKE) clean
	cd prdc; $(MAKE) clean
	cd r1d;  $(MAKE) clean
	cd rp;   $(MAKE) clean
	cd rpc;  $(MAKE) clean
	cd rpg;  $(MAKE) clean
	cd cpc;  $(MAKE) clean

# Clean unit tests
clean-tests:
	cd util/tests; $(MAKE) clean
	cd pscf/tests; $(MAKE) clean
	cd prdc/tests; $(MAKE) clean
	cd r1d/tests; $(MAKE) clean
	cd rp/tests/cpp;  $(MAKE) clean
	cd rp/tests/cuda;  $(MAKE) clean
	cd rpc/tests; $(MAKE) clean
	cd rpg/tests; $(MAKE) clean
	cd cpc/tests; $(MAKE) clean

# Remove all automatically generated files, recreate initial state
veryclean:
	cd util; $(MAKE) veryclean
	cd pscf; $(MAKE) veryclean
	cd prdc; $(MAKE) clean
	cd r1d; $(MAKE) clean
	cd rp;  $(MAKE) clean
	cd rpc; $(MAKE) clean
	cd rpg; $(MAKE) clean
	cd cpc; $(MAKE) clean
	rm -f util/config.mk
ifneq ($(BLD_DIR),$(SRC_DIR))
	rm -f util/makefile
	rm -f pscf/makefile
	rm -f prdc/makefile
	rm -f r1d/makefile
	rm -f rp/makefile
	rm -f rpc/makefile
	rm -f rpg/makefile
	rm -f cpc/makefile
	rm -f util/tests/makefile
	rm -f pscf/tests/makefile
	rm -f prdc/tests/makefile
	rm -f r1d/tests/makefile
	rm -f rp/tests/cpp/makefile
	rm -f rp/tests/cuda/makefile
	rm -f rpc/tests/makefile
	rm -f rpg/tests/makefile
	rm -f cpc/tests/makefile
	rm -f setopts
endif
	rm -f config.mk

# ======================================================================
