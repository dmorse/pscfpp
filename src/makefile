include config.mk

# ======================================================================
.PHONY: all fd1d pssp clean veryclean

all: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd fd1d; $(MAKE) all

fd1d: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd fd1d; $(MAKE) all

pssp: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd pssp; $(MAKE) all

# Run all unit tests
run-tests:
	cd util/tests; $(MAKE) all; $(MAKE) run
	cd pscf/tests; $(MAKE) all; $(MAKE) run
	cd fd1d/tests; $(MAKE) all; $(MAKE) run
	@cat util/tests/count > count
	@cat pscf/tests/count >> count
	@cat fd1d/tests/count >> count
	@echo " "
	@echo "Summary"
	@cat count
	@rm -f count

# Remove object (*.o), dependency (*.d) and library (*.a) files
clean:
	cd util; $(MAKE) clean
	cd pscf; $(MAKE) clean
	cd fd1d; $(MAKE) clean

# Clean unit tests
clean-tests:
	cd util/tests; $(MAKE) clean
	cd pscf/tests; $(MAKE) clean
	cd fd1d/tests; $(MAKE) clean

# Remove all automatically generated files, recreate initial state
veryclean:
	cd util; $(MAKE) veryclean
	cd pscf; $(MAKE) veryclean
	cd fd1d; $(MAKE) veryclean
	rm -f util/config.mk
	rm -f pscf/config.mk
	rm -f fd1d/config.mk
ifneq ($(BLD_DIR),$(SRC_DIR))
	rm -f util/makefile
	rm -f pscf/makefile
	rm -f fd1d/makefile
	rm -f util/tests/makefile
	rm -f pscf/tests/makefile
	rm -f fd1d/tests/makefile
	rm -f configure
endif
	rm -f config.mk

# ======================================================================
