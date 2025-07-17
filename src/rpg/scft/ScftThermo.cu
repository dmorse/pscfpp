/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ScftThermo.h"
#include <prdc/system/ScftReal.tpp>
#include <prdc/cuda/Reduce.h>
//#include <pscf/inter/Interaction.h>

namespace Pscf {
   namespace Prdc {
      template class ScftReal<1, Rpg::System<1> >;
      template class ScftReal<2, Rpg::System<2> >;
      template class ScftReal<3, Rpg::System<3> >;
   }
   namespace Rpg {

      /*
      * Constructor
      */
      template <int D>
      ScftThermo<D>::ScftThermo(System<D> const & system)
       : Base(system)
      {};

      /*
      * Inner product of r-grid fields.
      */
      template <int D>
      double ScftThermo<D>::innerProduct(FieldT const & A, 
                                         FieldT const & B) const
      {  return Cuda::Reduce::innerProduct(A, B); };

      // Explicit instantiation
      template class ScftThermo<1>;
      template class ScftThermo<2>;
      template class ScftThermo<3>;

   }
}
