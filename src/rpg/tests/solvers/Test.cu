#include <util/global.h>
#include "SolverTestComposite.h"

#include <test/TestRunner.h>
#include <test/CompositeTestRunner.h>
#include <util/param/BracketPolicy.h>

int main(int argc, char* argv[])
{
   BracketPolicy::set(BracketPolicy::Optional);

   SolverTestComposite runner;

   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
   }

   runner.run();
}
