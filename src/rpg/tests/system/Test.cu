#include <util/global.h>
#include "SystemTestComposite.h"

#include <util/param/BracketPolicy.h>
#include <test/TestRunner.h>
#include <test/CompositeTestRunner.h>

int main(int argc, char* argv[])
{
   BracketPolicy::set(BracketPolicy::Optional);
   SystemTestComposite runner;

   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
   }

   runner.run();
}
