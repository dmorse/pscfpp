/*
* This program runs all unit tests in the pspc/tests/cpu directory.
*/ 

#include "CpuTestComposite.h"
#include <util/global.h>

#include <test/TestRunner.h>
#include <test/CompositeTestRunner.h>

int main(int argc, char* argv[])
{
   CpuTestComposite runner;

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
