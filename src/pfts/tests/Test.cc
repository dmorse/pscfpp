/*
* This program runs all unit tests in the pfts/tests directory.
*/ 

#include <util/global.h>
#include "PftsTestComposite.h"

#include <test/CompositeTestRunner.h>

using namespace Pfts;
using namespace Util;

int main(int argc, char* argv[])
{
   PftsTestComposite runner;

   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
    }
   runner.run();
}
