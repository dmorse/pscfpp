include src/config.mk
# ==========================================================================
.PHONY: all fd1d test clean veryclean \
        html clean-html

# ==============================================================================
# Main build targets

all:
	cd bld; $(MAKE) all

fd1d:
	cd bld; $(MAKE) fd1d

# ==========================================================================
# Test targets

test:
	@cd src/util/tests; $(MAKE) all; $(MAKE) run
	@cd src/pscf/tests; $(MAKE) all; $(MAKE) run
	@cd src/fd1d/tests; $(MAKE) all; $(MAKE) run
	@cat src/util/tests/count >> count
	@cat src/chem/tests/count >> count
	@echo " "
	@echo "Summary"
	@cat count
	@rm -f count

# =========================================================================
# Clean targets

clean:
	cd bld; $(MAKE) clean
	cd src; $(MAKE) clean

clean-tests:
	cd src/; $(MAKE) clean-tests

clean-bin:
	-rm -f $(BIN_DIR)/makeDep
	-rm -f $(BIN_DIR)/pscf*
 
veryclean:
	make clean-bin
	cd bld; $(MAKE) veryclean
	rm bld/makefile
	cd src; $(MAKE) veryclean
	cd doc; $(MAKE) clean

# =========================================================================
# HTML Documentation
 
html:
	cd doc; $(MAKE) html

clean-html:
	cd doc; $(MAKE) clean

# ==========================================================================
