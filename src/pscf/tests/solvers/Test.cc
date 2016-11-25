/*
* This program runs all unit tests in the pscf/tests/solvers directory.
*/ 

#include <util/global.h>
#include "SolversTestComposite.h"

#include <test/CompositeTestRunner.h>

using namespace Pscf;
using namespace Util;

int main(int argc, char* argv[])
{
   SolversTestComposite runner;

   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
    }
   runner.run();
}
