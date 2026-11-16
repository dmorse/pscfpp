/*
* This program runs all unit tests in the rpc/tests/field directory.
*/ 

//#include "RFieldTest.h"
#include "FieldTestComposite.h"
#include <util/param/BracketPolicy.h>
//#include <util/global.h>

#include <test/TestRunner.h>
#include <test/CompositeTestRunner.h>

int main(int argc, char* argv[])
{
   FieldTestComposite runner;
   BracketPolicy::set(BracketPolicy::Optional);

   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
   }

   runner.run();
}
