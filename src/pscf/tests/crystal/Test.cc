/*
* This program runs all unit tests in the pscf/tests/math directory.
*/ 

#include <util/global.h>
#include "UnitCellTest.h"

#include <test/TestRunner.h>
#include <test/CompositeTestRunner.h>

using namespace Pscf;
using namespace Util;

int main(int argc, char* argv[])
{
   TEST_RUNNER(UnitCellTest) runner;

   // MathTestComposite runner;

   #if 0
   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
   }
   #endif
   runner.run();
}
