#ifndef RPC_CONCENTRATION_DERIVATIVE_H
#define RPC_CONCENTRATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"
#include <rp/fts/analyzer/ConcentrationDerivative.h>
#include <rpc/system/Types.h>

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Evaluate the derivative of H with respect to chi.
   *
   * Specializations of this template are derived from specializations of 
   * the base class template Rp::ConcentrationDerivative, and inherit their 
   * entire public interface and almost all of their source code from this 
   * base class. 
   *
   * \see Rp::ConcentrationDerivative
   * \see \ref rp_ConcentrationDerivative_page "Manual Page"
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class ConcentrationDerivative
    : public Rp::ConcentrationDerivative< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      ConcentrationDerivative(Simulator<D>& simulator, Rp::System<D, Rpc::Types<D> >& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class ConcentrationDerivative<1, Rpc::Types<1> >;
      extern template class ConcentrationDerivative<2, Rpc::Types<2> >;
      extern template class ConcentrationDerivative<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class ConcentrationDerivative<1>;
      extern template class ConcentrationDerivative<2>;
      extern template class ConcentrationDerivative<3>;
   }
}
#endif
