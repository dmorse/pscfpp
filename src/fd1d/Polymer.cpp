/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"

namespace Pscf {
namespace Fd1d 
{ 

   Polymer::Polymer()
   {  setClassName("Polymer"); }

   Polymer::~Polymer()
   {}

   void Polymer::setPhi(double phi)
   {
      UTIL_CHECK(phi >= 0.0);  
      UTIL_CHECK(phi <= 1.0);  
      phi_ = phi; 
   }

}
}
