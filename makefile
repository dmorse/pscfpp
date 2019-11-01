include src/config.mk
# =========================================================================
.PHONY:  html clean clean-tests clean-bin clean-html veryclean 

# =========================================================================
# HTML Documentation
 
html:
	cd doc; $(MAKE) html

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
	cd doc; $(MAKE) clean

veryclean:
	make clean-bin
	rm -f $(BIN_DIR)/makeDep*
	cd bld; $(MAKE) veryclean
	rm bld/makefile
	cd src; $(MAKE) veryclean
	cd doc; $(MAKE) clean



# ==========================================================================
