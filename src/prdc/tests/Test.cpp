/*
* This program runs all unit tests in the prdc/tests directory.
*/ 

#include <test/CompositeTestRunner.h>

#include "crystal/CrystalTestComposite.h"

#include <util/param/BracketPolicy.h>
#include <util/global.h>

TEST_COMPOSITE_BEGIN(PscfNsTestComposite)
addChild(new CrystalTestComposite, "crystal/");
TEST_COMPOSITE_END

using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Util;

int main(int argc, char* argv[])
{

   BracketPolicy::set(BracketPolicy::Optional);

   try {
   
      PscfNsTestComposite runner;

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

      std::cerr << "Uncaught exception in prdc/tests/Test.cc" << std::endl;
      return 1;

   }
}
