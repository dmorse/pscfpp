#ifndef RPC_PERTURBATION_DERIVATIVE_H
#define RPC_PERTURBATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"
#include <rp/fts/analyzer/PerturbationDerivative.h>
#include <rpc/system/Types.h>

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Evaluate derivative of H w/ respect to perturbation parameter lambda.
   *
   * Specializations of this template are derived from specializations of 
   * the base class template Rp::PerturbationDerivative, and inherit their 
   * entire public interface and almost all of their source code from this 
   * base class. See the documentation for this base class template for
   * details. 
   *
   * \see Rp::PerturbationDerivative
   * \see rp_PerturbationDerivative_page "Manual Page"
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class PerturbationDerivative 
     : public Rp::PerturbationDerivative< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      PerturbationDerivative(Simulator<D>& simulator, Rp::System<D, Rpc::Types<D> >& system);

   };

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class PerturbationDerivative< 1, Rpc::Types<1> >;
      extern template class PerturbationDerivative< 2, Rpc::Types<2> >;
      extern template class PerturbationDerivative< 3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class PerturbationDerivative<1>;
      extern template class PerturbationDerivative<2>;
      extern template class PerturbationDerivative<3>;
   }
}
#endif
