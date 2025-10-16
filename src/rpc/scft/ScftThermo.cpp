/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ScftThermo.h"
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <prdc/system/ScftThermoTmpl.tpp>
#include <pscf/inter/Interaction.h>

namespace Pscf {

   namespace Prdc {
      template class ScftThermoTmpl<1, Rpc::System<1> >;
      template class ScftThermoTmpl<2, Rpc::System<2> >;
      template class ScftThermoTmpl<3, Rpc::System<3> >;
   }

   namespace Rpc {

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
      double ScftThermo<D>::innerProduct(RFieldT const & A, 
                                         RFieldT const & B) const
      {
         const int meshSize = Base::domain().mesh().size();
         UTIL_CHECK(meshSize == A.capacity())
         UTIL_CHECK(meshSize == B.capacity())
         double sum = 0.0;
         for (int k = 0; k < meshSize; ++k) {
            sum += A[k]*B[k];
         }
         return sum;
      };

      // Explicit instantiation
      template class ScftThermo<1>;
      template class ScftThermo<2>;
      template class ScftThermo<3>;
   }
}
