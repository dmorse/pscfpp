/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CUT.h>
#include <rp/fts/perturbation/Perturbation.tpp> // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Perturbation<1,CUT>;
      template class Perturbation<2,CUT>;
      template class Perturbation<3,CUT>;
   }
}
