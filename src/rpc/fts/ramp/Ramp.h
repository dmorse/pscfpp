#ifndef RPC_RAMP_H
#define RPC_RAMP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/ramp/Ramp.h>            // base class template
#include <rpc/system/Types.h>            // base class template argument

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /**
   * Ramp that varies parameters linearly with index.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::Ramp, and inherit their 
   * entire public interface and almost all of their source code from this
   * base class.  See the documentation of this base class template for 
   * details. 
   *
   * \see \ref Rp::Ramp
   * \see \ref psfts_ramp_page
   * \ingroup Rpc_Fts_Ramp_Module
   */
   template <int D>
   class Ramp : public Rp::Ramp<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator
      */
      Ramp(Rp::Simulator<D, Rpc::Types<D> >& simulator);

      /**
      * Destructor.
      */
      virtual ~Ramp() = default;

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Ramp<1, Rpc::Types<1> >;
      extern template class Ramp<2, Rpc::Types<2> >;
      extern template class Ramp<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class Ramp<1>;
      extern template class Ramp<2>;
      extern template class Ramp<3>;
   }
}
#endif
