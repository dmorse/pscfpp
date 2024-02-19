/*
* This program runs all unit tests in the r1d/tests directory.
*/ 

#include <util/global.h>
#include <util/param/BracketPolicy.h>

#include "R1dTestComposite.h"

#include <test/CompositeTestRunner.h>

using namespace R1d;
using namespace Util;

int main(int argc, char* argv[])
{

   BracketPolicy::set(BracketPolicy::Optional);

   try {

      R1dTestComposite runner;
  
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

      std::cerr << "Uncaught exception in r1d/tests/Test.cc" << std::endl;
      return 1;

   }
}
