/*
* This program runs all unit tests in the fd1d/tests directory.
*/ 

#include <util/global.h>
#include "Fd1dTestComposite.h"

#include <test/CompositeTestRunner.h>

using namespace Fd1d;
using namespace Util;

int main(int argc, char* argv[])
{
   Fd1dTestComposite runner;

   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
   }
   runner.run();
}
