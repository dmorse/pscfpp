#ifndef RPC_CHI_DERIVATIVE_H
#define RPC_CHI_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"                // indirect base class
#include <rp/fts/analyzer/BinaryChiDerivative.h>  // base class template
#include <rpc/system/Types.h>               // base template argument

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /**
   * Evaluate the derivative of H with respect to chi.
   *
   * Specializations of this class template are derived from specializations 
   * of the base class template Rp::BinaryChiDerivative, and inherit their public 
   * interface and almost all of their source code from this base class.
   *
   * \see Rp::BinaryChiDerivative
   * \see \ref rp_BinaryChiDerivative_page "Manual Page"
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class BinaryChiDerivative : public Rp::BinaryChiDerivative< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      BinaryChiDerivative(Rp::Simulator<D, Rpc::Types<D> >& simulator, Rp::System<D, Rpc::Types<D> >& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BinaryChiDerivative<1, Rpc::Types<1> >;
      extern template class BinaryChiDerivative<2, Rpc::Types<2> >;
      extern template class BinaryChiDerivative<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class BinaryChiDerivative<1>;
      extern template class BinaryChiDerivative<2>;
      extern template class BinaryChiDerivative<3>;
   }
}
#endif
