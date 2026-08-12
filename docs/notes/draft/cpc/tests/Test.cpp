/*
* This program runs all unit tests in the pscf/tests directory.
*/ 

#include <test/CompositeTestRunner.h>

#include "solvers/SolverTestComposite.h"
#include "field/FieldTestComposite.h"

#include <util/param/BracketPolicy.h>
#include <util/global.h>

TEST_COMPOSITE_BEGIN(CpcNsTestComposite)
addChild(new SolverTestComposite, "solvers/");
addChild(new FieldTestComposite, "field/");
TEST_COMPOSITE_END

using namespace Util;
using namespace Pscf;

int main(int argc, char* argv[])
{

   BracketPolicy::set(BracketPolicy::Optional);

   try {

      CpcNsTestComposite runner;

      // Add any file prefix given as command line argument
      if (argc > 2) {
         UTIL_THROW("Too many arguments");
      }
      if (argc == 2) {
         runner.addFilePrefix(argv[1]);
      }
   
      // Run all unit test methods
      int failures = runner.run();

      if (failures != 0) {
         failures = 1;
      }
      return failures;

   } catch (...) {

      std::cerr << "Uncaught exception in cpc/tests/Test.cc" << std::endl;
      return 1;

   }
}
