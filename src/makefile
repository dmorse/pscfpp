# ======================================================================
.PHONY: all clean veryclean

all: 
	./configure -m0
	cd util; $(MAKE) all
	cd chem; $(MAKE) all

# Remove object (*.o), dependency (*.d) and library (*.a) files
clean:
	cd util; $(MAKE) clean
	cd chem; $(MAKE) clean

# Clean unit tests
clean-tests:
	cd util/tests; $(MAKE) clean
	cd chem/tests; $(MAKE) clean

# Remove all automatically generated files, recreate initial state
veryclean:
	cd util; $(MAKE) veryclean
	cd chem; $(MAKE) veryclean
	rm -f util/config.mk
	rm -f chem/config.mk
	rm -f config.mk

# ======================================================================
