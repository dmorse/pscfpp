/*
* This program runs all unit tests in the chem/tests directory.
*/ 

#include <util/global.h>
#include "ChemTestComposite.h"

#include <test/CompositeTestRunner.h>

using namespace Pfts;
using namespace Util;

int main(int argc, char* argv[])
{
   ChemTestComposite runner;

   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
    }
   runner.run();
}
