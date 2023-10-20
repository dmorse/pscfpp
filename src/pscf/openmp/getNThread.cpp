/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "getNThread.h"
#include <util/global.h>

#include <unistd.h>
#include <cstring>

using namespace Util;

namespace Pscf { 
namespace Prdc { 

   /*
   * Extract integer argument of -d option.
   */
   int getNThread(int argc, char **argv)
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
                  if (option[1] == 't') {
                     found = true;
                  }
               }
               if (found && length > 2) {
                  // Option of form -d# with no space
                  arg = &option[2];
                  done = true;
               }
            }
         } else {
            arg = option;
            done = true;
         }
         ++i;
      }
      int nThread = 0;
      if (done) {
         nThread = atoi(arg); 
      }
      return nThread; 
   }

}
}
