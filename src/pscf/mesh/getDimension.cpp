/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "getDimension.h"
#include <util/global.h>

//#include <iostream>
//#include <string.h>
#include <unistd.h>

using namespace Util;

namespace Pscf { 

   /*
   * Extract integer argument of -d option.
   */
   int getDimension(int argc, char **argv)
   {
      // No options beyond program name
      if (argc < 2) return 0;

      int length;
      char* option = 0;
      char* arg = 0;

      bool next = false;
      bool done = false;
      int i = 1;
      while (!done && i < argc) {
         option = argv[i];
         length = strlen(option);
         if (!next) {
            if (length > 1) {
               if (option[0] == '-') {
                   if (option[1] == 'd') {
                      next = true;
                   }
               }
               if (next && length > 2) {
                  if (length == 3) {
                     // Option of form -d# with no space
                     arg = &option[2];
                     next = false;
                     done = true;
                  } else {
                     // Option is too long
                     return 0;
                  }
               }
            }
         } else {
            if (length == 1) {
               // Single character option
               arg = option;
               next = false;
               done = true;
            } else {
               // Option is too long
               return 0;
            }
         }
         ++i;
      }
      UTIL_CHECK(1 == strlen(arg));

      // Convert arg string to integer D
      // int D;
      // std::sscanf(arg, "%d", &D);
      int D = atoi(arg);
      UTIL_CHECK(D > 0 && D < 4);

      return D;
   }

}
