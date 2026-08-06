/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CUT.h>    // backend type class
#include <rp/field/CFields.tpp>   // clas template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class CFields<1,CUT>;
      template class CFields<2,CUT>;
      template class CFields<3,CUT>;
   } 
}
