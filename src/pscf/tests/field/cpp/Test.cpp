/*
* This program runs all unit tests in the pscf/tests/field/cpp directory.
*/ 

#include "FieldCppTestComposite.h"
#include <util/global.h>

#include <test/TestRunner.h>
#include <test/CompositeTestRunner.h>

int main(int argc, char* argv[])
{
   FieldCppTestComposite runner;

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
