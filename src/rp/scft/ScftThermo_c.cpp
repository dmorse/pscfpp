/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cpu/Reduce.h>

#include <pscf/backends/CPT.h>
#include <rp/scft/ScftThermo.tpp>     // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ScftThermo<1,CPT>;
      template class ScftThermo<2,CPT>;
      template class ScftThermo<3,CPT>;
   }
}
