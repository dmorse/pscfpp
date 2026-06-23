#ifndef RPC_LINEAR_RAMP_H
#define RPC_LINEAR_RAMP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/ramp/LinearRamp.h>      // direct base class template
#include <rpc/system/Types.h>            // base class template argument
#include <rpc/fts/ramp/RampParameter.h>  // base class member
#include <rpc/fts/ramp/Ramp.h>           // indirect base class

namespace Pscf {
namespace Rpc {

   // Forward declaration
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Ramp that varies parameters linearly with index.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::LinearRamp, and
   * inherit their public interface and almost all of their source code
   * from this base class.  See the documentation of this base class 
   * template for details. 
   *
   * \see \ref Rp::LinearRamp
   * \see \ref rp_LinearRamp_page "Manual Page"
   * \ingroup Rpc_Fts_Ramp_Module
   */
   template <int D>
   class LinearRamp : public Rp::LinearRamp<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator
      */
      LinearRamp(Rp::Simulator<D, Rpc::Types<D> >& simulator);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class LinearRamp<1, Rpc::Types<1> >;
      extern template class LinearRamp<2, Rpc::Types<2> >;
      extern template class LinearRamp<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class LinearRamp<1>;
      extern template class LinearRamp<2>;
      extern template class LinearRamp<3>;
   }
}
#endif
