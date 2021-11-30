include src/config.mk
# =========================================================================
.PHONY:  html clean clean-tests clean-bin clean-html veryclean 

# =========================================================================
# Main targets 
 
all-cpu:
	cd bld; $(MAKE) all-cpu

fd1d:
	cd bld; $(MAKE) fd1d

pspc:
	cd bld; $(MAKE) pspc


# =========================================================================
# HTML Documentation
 
html:
	cd docs; $(MAKE) html

# =========================================================================
# Clean targets

clean:
	cd bld; $(MAKE) clean
	cd src; $(MAKE) clean

clean-tests:
	cd src/; $(MAKE) clean-tests
	cd bld/; $(MAKE) clean-tests

clean-bin:
	rm -f $(BIN_DIR)/pscf*
 
clean-html:
	cd docs; $(MAKE) clean

veryclean:
	make clean-bin
	rm -f $(BIN_DIR)/makeDep*
	cd bld; $(MAKE) veryclean
	rm bld/makefile
	cd src; $(MAKE) veryclean
	cd docs; $(MAKE) clean

# ==========================================================================
