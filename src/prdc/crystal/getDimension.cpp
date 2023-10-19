/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "getDimension.h"
#include <util/global.h>

#include <unistd.h>
#include <cstring>

using namespace Util;

namespace Pscf { 
namespace Prdc { 

   /*
   * Extract integer argument of -d option.
   */
   int getDimension(int argc, char **argv)
   {
      // No arguments other than program name
      if (argc < 2) {
         UTIL_THROW("No command options found");
      }

      int length;
      char* option = 0;
      char* arg = 0;

      bool found = false;  // Has the -d option been found
      bool done = false;   // Has a valid parameter been found
      int i = 1;
      while (!done && i < argc) {
         option = argv[i];
         length = strlen(option);
         if (!found) {
            if (length > 1) {
               if (option[0] == '-') {
                   if (option[1] == 'd') {
                      found = true;
                   }
               }
               if (found && length > 2) {
                  if (length == 3) {
                     // Option of form -d# with no space
                     arg = &option[2];
                     done = true;
                  } else {
                     // Option is too long
                     std::cout 
                         << "Invalid parameter of command option -d:"
                         << " |" << option << "|" << std::endl;
                     UTIL_THROW("Invalid parameter of command option -d");
                  }
               }
            }
         } else {
            if (length == 1) {
               // Single character option parameter
               arg = option;
               done = true;
            } else {
               // Option is too long
               std::cout << "Invalid parameter of command option -d:"
                         << " |" << option << "|" << std::endl;
               UTIL_THROW("Invalid parameter of command option -d");
            }
         }
         ++i;
      }

      // Check status
      if (!done) {
         if (found) {
            UTIL_THROW("Missing parameter for option -d");
         } else {
            UTIL_THROW("Missing required option -d");
         }
      }
      UTIL_CHECK(1 == strlen(arg));

      // Convert arg string to integer D
      int D = atoi(arg);
      if (D > 0 && D < 4) {
         return D;
      } else {
         std::cout << "Invalid parameter of command line option -d:"
                   << " |" << arg << "|" << std::endl;
         std::cout << "Value of attempted conversion to integer: " 
                   << D << std::endl;
         UTIL_THROW("Invalid parameter of command line option -d");
      }

   }

}
}
