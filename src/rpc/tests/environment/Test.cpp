/*
* This program runs all unit tests in the rpc/tests/environment directory.
*/ 

#include "EnvironmentTestComposite.h"

#include <test/TestRunner.h>
#include <test/CompositeTestRunner.h>

int main(int argc, char* argv[])
{
   EnvironmentTestComposite runner;

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
