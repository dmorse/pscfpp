/*
* This program runs all unit tests in the pscf/tests directory.
*/ 

#include <test/CompositeTestRunner.h>

#include "field/FieldTestComposite.h"
#include "solvers/SolverTestComposite.h"
#include "iterator/IteratorTest.h"
#include <util/global.h>

TEST_COMPOSITE_BEGIN(PsspNsTestComposite)
addChild(new FieldTestComposite, "field/");
addChild(new SolverTestComposite, "solvers/");
//addChild(new TEST_RUNNER(IteratorTest), "iterator/");
TEST_COMPOSITE_END

using namespace Pscf;
using namespace Util;

int main(int argc, char* argv[])
{
   PsspNsTestComposite runner;

   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
    }
   runner.run();
}
