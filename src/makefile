include config.mk

.PHONY: all all-cpu util pscf fd1d pspc pspg \
        test test-cpu test-gpu \
        clean clean-tests veryclean
# ======================================================================
# Main build targets

# Build all CPU and GPU programs
all:
	$(MAKE) all-cpu
ifdef PSCF_CUDA
	$(MAKE) pspg
endif

# Build all programs that uses exclusively CPU hardware
all-cpu: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd fd1d; $(MAKE) all
	cd prdc; $(MAKE) all
	cd pspc; $(MAKE) all

# Build code in Util names (general scientific utilities)
util: 
	cd util; $(MAKE) all

# Build common code for PSCF in Util and Pscf namespaces
pscf: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all

# Build common code for PSCF in Util, Pscf, and Prdc namespaces
prdc: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd prdc; $(MAKE) all

# Build pscf_fd 1D finite difference SCFT program (install in BIN_DIR)
pscf_fd: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd fd1d; $(MAKE) all

# Build pscf_pc CPU code for periodic structures (install in BIN_DIR)
pscf_pc: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd prdc; $(MAKE) all
	cd pspc; $(MAKE) all

ifdef PSCF_CUDA
# Build pscf_pc GPU code for periodic structures (install in BIN_DIR)
pscf_pg:
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd prdc; $(MAKE) all
	cd pspg; $(MAKE) all
endif

# ======================================================================
# Test targets

# Run all CPU and GPU unit tests
test:
	@cd util/tests; $(MAKE) all; $(MAKE) run-quiet
	@cd pscf/tests; $(MAKE) all; $(MAKE) run-quiet
	@cd prdc/tests; $(MAKE) all; $(MAKE) run-quiet
	@cd fd1d/tests; $(MAKE) all; $(MAKE) run-quiet
	@cd pspc/tests; $(MAKE) all; $(MAKE) run-quiet
ifdef PSCF_CUDA
	@cd pspg/tests; $(MAKE) all; $(MAKE) run-quiet
endif
	@cat util/tests/count > count
	@cat pscf/tests/count >> count
	@cat prdc/tests/count >> count
	@cat fd1d/tests/count >> count
	@cat pspc/tests/count >> count
ifdef PSCF_CUDA
	@cat pspg/tests/count >> count
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
	@cd fd1d/tests; $(MAKE) all; $(MAKE) run-quiet
	@cd pspc/tests; $(MAKE) all; $(MAKE) run-quiet
	@cat util/tests/count > count
	@cat pscf/tests/count >> count
	@cat prdc/tests/count >> count
	@cat fd1d/tests/count >> count
	@cat pspc/tests/count >> count
	@echo " "
	@echo "Summary"
	@cat count
	@rm -f count

test-gpu:
	@cd pspg/tests; @$(MAKE) all; $(MAKE) run-quiet

# ======================================================================
# Clean targets

# Remove object (*.o), dependency (*.d) and library (*.a) files
clean:
	cd util; $(MAKE) clean
	cd pscf; $(MAKE) clean
	cd fd1d; $(MAKE) clean
	cd pspc; $(MAKE) clean
	cd pspg; $(MAKE) clean

# Clean unit tests
clean-tests:
	cd util/tests; $(MAKE) clean
	cd pscf/tests; $(MAKE) clean
	cd fd1d/tests; $(MAKE) clean
	cd pspc/tests; $(MAKE) clean
	cd pspg/tests; $(MAKE) clean

# Remove all automatically generated files, recreate initial state
veryclean:
	cd util; $(MAKE) veryclean
	cd pscf; $(MAKE) veryclean
	cd prdc; $(MAKE) veryclean
	cd fd1d; $(MAKE) veryclean
	cd pspc; $(MAKE) veryclean
	cd pspg; $(MAKE) veryclean
	rm -f util/config.mk
	rm -f pscf/config.mk
	rm -f prdc/config.mk
	rm -f fd1d/config.mk
	rm -f pspc/config.mk
	rm -f pspg/config.mk
ifneq ($(BLD_DIR),$(SRC_DIR))
	rm -f util/makefile
	rm -f pscf/makefile
	rm -f prdc/makefile
	rm -f fd1d/makefile
	rm -f pspc/makefile
	rm -f pspg/makefile
	rm -f util/tests/makefile
	rm -f pscf/tests/makefile
	rm -f prdc/tests/makefile
	rm -f fd1d/tests/makefile
	rm -f pspc/tests/makefile
	rm -f setopts
endif
	rm -f config.mk

# ======================================================================
