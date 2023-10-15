/*
* This program runs all unit tests in the prdc/tests directory.
*/ 

#include <test/CompositeTestRunner.h>

#include "crystal/CrystalTestComposite.h"
#include "cpu/CpuTestComposite.h"
#ifdef PSCF_CUDA
#include "cuda/CudaTestComposite.h"
#endif

#include <util/param/BracketPolicy.h>
#include <util/global.h>

TEST_COMPOSITE_BEGIN(PrdcTestComposite)
addChild(new CrystalTestComposite, "crystal/");
addChild(new CpuTestComposite, "cpu/");
#ifdef PSCF_CUDA
addChild(new CudaTestComposite, "cuda/");
#endif
TEST_COMPOSITE_END

using namespace Util;
using namespace Pscf;

int main(int argc, char* argv[])
{

   BracketPolicy::set(BracketPolicy::Optional);

   try {
   
      PrdcTestComposite runner;

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
