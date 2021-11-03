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
   try {

      Fd1dTestComposite runner;
  
      // Optionally add file prefix given as command line argument 
      if (argc > 2) {
         UTIL_THROW("Too many arguments");
      }
      if (argc == 2) {
         runner.addFilePrefix(argv[1]);
      }
   
      // Run all unit test functions
      int failures = runner.run();

      return (failures != 0);

   } catch (...) {

      std::cerr << "Uncaught exception in fd1d/tests/Test.cc" << std::endl;
      return 1;

   }
}
