/*
* This program runs all unit tests in the pfts/tests directory.
*/ 

#include <util/global.h>
#include "chem/ChemTestComposite.h"

#include <test/CompositeTestRunner.h>

using namespace Pfts;
using namespace Util;

TEST_COMPOSITE_BEGIN(PftsNsTestComposite)
addChild(new ChemTestComposite, "chem/");
TEST_COMPOSITE_END

int main(int argc, char* argv[])
{
   PftsNsTestComposite runner;

   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
    }
   runner.run();
}
