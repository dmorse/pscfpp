/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CPT.h> 
#include <rp/fts/ramp/LinearRamp.tpp> 


// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LinearRamp<1,CPT>;
      template class LinearRamp<2,CPT>;
      template class LinearRamp<3,CPT>;
   }
}
