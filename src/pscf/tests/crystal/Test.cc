/*
* This program runs all unit tests in the pscf/tests/math directory.
*/ 

#include <util/global.h>
#include "CrystalTestComposite.h"

#include <test/TestRunner.h>
#include <test/CompositeTestRunner.h>

using namespace Pscf;
using namespace Util;

int main(int argc, char* argv[])
{
   #if 0
   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
   }
   #endif

   CrystalTestComposite runner;
   runner.run();

   #if 0
   TEST_RUNNER(UnitCellTest) runner1;
   runner1.run();

   TEST_RUNNER(SpaceSymmetryTest) runner2;
   runner2.run();
   #endif

}
