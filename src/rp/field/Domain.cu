/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

//#include <rp/field/FieldIo_u.h> 
#include <rp/field/Domain.tpp> 

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Domain<1,CUT>;
      template class Domain<2,CUT>;
      template class Domain<3,CUT>;
   }
}
