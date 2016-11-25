/*
* This program runs all unit tests in the pscf/tests directory.
*/ 

#include <test/CompositeTestRunner.h>

#include "math/MathTestComposite.h"
#include "chem/ChemTestComposite.h"
#include "solvers/SolversTestComposite.h"
#include "inter/InterTestComposite.h"
#include <util/global.h>

TEST_COMPOSITE_BEGIN(PscfNsTestComposite)
addChild(new MathTestComposite, "math/");
addChild(new ChemTestComposite, "chem/");
addChild(new SolversTestComposite, "solvers/");
addChild(new InterTestComposite, "inter/");
TEST_COMPOSITE_END

using namespace Pscf;
using namespace Util;

int main(int argc, char* argv[])
{
   PscfNsTestComposite runner;

   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
    }
   runner.run();
}
