# ======================================================================
.PHONY: all clean veryclean

all: 
	cd util; $(MAKE) all
	cd pfts; $(MAKE) all

# Remove object (*.o), dependency (*.d) and library (*.a) files
clean:
	cd util; $(MAKE) clean
	cd pfts; $(MAKE) clean

# Clean unit tests
clean-tests:
	cd util/tests; $(MAKE) clean
	cd pfts/tests; $(MAKE) clean

# Remove all automatically generated files, recreate initial state
veryclean:
	cd util; $(MAKE) veryclean
	cd pfts; $(MAKE) veryclean
	rm -f util/config.mk
	rm -f pfts/config.mk
	rm -f config.mk

# ======================================================================
