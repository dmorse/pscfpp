/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/system/System.h>
#include <rp/field/Domain.h>

#include <rp/system/SystemConstRef.tpp>

// Explicit initialization definitions
namespace Pscf {
   namespace Rp {
      template class SystemConstRef<1,CPT>;
      template class SystemConstRef<2,CPT>;
      template class SystemConstRef<3,CPT>;
   }
} 
