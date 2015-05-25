include src/config.mk
# ==========================================================================
.PHONY: all test clean veryclean \
        html clean-html

# ==============================================================================
# Main build targets

all:
	cd src; $(MAKE) all

# ==========================================================================
# Test targets

test:
	@cd src/util/tests; $(MAKE) all; $(MAKE) run
	@cd src/chem/tests; $(MAKE) all; $(MAKE) run
	@cat src/util/tests/count >> count
	@cat src/chem/tests/count >> count
	@echo " "
	@echo "Summary"
	@cat count
	@rm -f count

# =========================================================================
# Clean targets

clean:
	cd src; $(MAKE) clean

clean-tests:
	cd src/; $(MAKE) clean-tests

clean-bin:
	-rm -f $(BIN_DIR)/makeDep
	-rm -f $(BIN_DIR)/pft*
 
veryclean:
	make clean-bin
	cd src; $(MAKE) veryclean
#	cd doc; $(MAKE) clean

# =========================================================================
# HTML Documentation
 
#html:
#	cd doc; $(MAKE) html

#clean-html:
#	cd doc; $(MAKE) clean

# ==========================================================================
