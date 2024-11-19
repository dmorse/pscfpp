/*
* This program runs all unit tests in the rpg/tests/fts directory.
*/ 

#include <util/global.h>
#include "SimulatorTest.h"

#include <test/TestRunner.h>
#include <test/CompositeTestRunner.h>

int main(int argc, char* argv[])
{
   TEST_RUNNER(SimulatorTest) runner;

   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
   }

   runner.run();
}
