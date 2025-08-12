/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "NanException.h"
#include <sstream>

namespace Pscf
{

   using namespace Util;

   /*
   * Constructor, with function name.
   */
   NanException::NanException(const char *function, const char *file, 
                              int line, int echo) 
    : Exception(function, "required numerical parameter has value of NaN",
                file, line, echo)
   {}

   /*
   * Constructor, with no function name.
   */
   NanException::NanException(const char *file, int line, int echo) 
    : Exception("required numerical parameter has value of NaN",
                file, line, echo)
   {}

   /*
   * Destructor.
   */
   NanException::~NanException()
   {}

}