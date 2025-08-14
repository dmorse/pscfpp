/*
* This program runs all unit tests in the pscf/tests directory.
*/ 

#include <test/CompositeTestRunner.h>

#include "math/MathTestComposite.h"
#include "chem/ChemTestComposite.h"
#include "solvers/SolversTestComposite.h"
#include "inter/InterTestComposite.h"
#include "floryHuggins/FloryHugginsTestComposite.h"
#include "mesh/MeshTestComposite.h"
#ifdef PSCF_CUDA
#include "cuda/CudaTestComposite.h"
#endif

#include <util/param/BracketPolicy.h>
#include <util/global.h>

TEST_COMPOSITE_BEGIN(PscfNsTestComposite)
addChild(new MathTestComposite, "math/");
addChild(new ChemTestComposite, "chem/");
addChild(new SolversTestComposite, "solvers/");
addChild(new InterTestComposite, "inter/");
addChild(new FloryHugginsTestComposite, "floryHuggins/");
addChild(new MeshTestComposite, "mesh/");
#ifdef PSCF_CUDA
addChild(new CudaTestComposite, "cuda/");
#endif
TEST_COMPOSITE_END

using namespace Pscf;
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

      std::cerr << "Uncaught exception in pscf/tests/Test.cc" << std::endl;
      return 1;

   }
}
