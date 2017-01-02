SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk

#-----------------------------------------------------------------------
EXES=example1 example2

ifdef UTIL_MPI
EXES+=example3
endif

#-----------------------------------------------------------------------
all: $(EXES)

clean:
	rm -f $(EXES)
#-----------------------------------------------------------------------
# Path(s) to search for header files. 
TEST_DIR=$(SRC_DIR)/test
INCLUDES= -I$(TEST_DIR)

# Extra dependencies for all source files
ALLDEPS= $(SRC_DIR)/config.mk \
        $(TEST_DIR)/UnitTest.h $(TEST_DIR)/TestRunner.h \
        $(TEST_DIR)/UnitTestRunner.h $(TEST_DIR)/CompositeTestRunner.h
#-----------------------------------------------------------------------

example1: example1.cpp $(ALLDEPS)
	$(CXX) $(CPPFLAGS) $(TESTFLAGS) $(INCLUDES) -o $@ $<
	
example2: example2.cpp $(ALLDEPS)
	$(CXX) $(CPPFLAGS) $(TESTFLAGS) $(INCLUDES) -o $@ $<

example3: example3.cpp $(ALLDEPS)
	$(CXX) $(CPPFLAGS) $(TESTFLAGS) $(INCLUDES) -o $@ $<
	
