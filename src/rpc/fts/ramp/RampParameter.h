#ifndef RPC_RAMP_PARAMETER_H
#define RPC_RAMP_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/ramp/RampParameter.h>   // direct base class template
#include <rpc/system/Types.h>            // base class template argument
#include <rpc/fts/ramp/Ramp.h>           // indirect base class

namespace Pscf {
namespace Rpc {

   template <int D> class Simulator;

   using namespace Util;

   /**
   * Class for storing data about an individual ramp parameter.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::RampParameter, and
   * inherit their public interface and almost all of their source code
   * from this base class.  See the documentation of this base class 
   * template for details. 
   *
   * \see \ref Rp::RampParameter
   * \see \ref rp_LinearRamp_page "Manual Page"
   * \ingroup Rpc_Fts_Ramp_Module
   */
   template <int D>
   class RampParameter : public Rp::RampParameter<D, Types<D> >
   {

   public:

      /**
      * Default constructor.
      */
      RampParameter();

      /**
      * Constructor that stores a pointer to parent Simulator.
      *
      * \param simulator  parent Simulator
      */
      RampParameter(Simulator<D>& simulator);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class RampParameter<1, Rpc::Types<1> >;
      extern template class RampParameter<2, Rpc::Types<2> >;
      extern template class RampParameter<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class RampParameter<1>;
      extern template class RampParameter<2>;
      extern template class RampParameter<3>;
   }
}
#endif
