#ifndef RPC_SCFT_THERMO_H
#define RPC_SCFT_THERMO_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/ScftThermo.h>    // base class template
#include <rpc/system/System.h>     // template argument

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * Computes SCFT free energies.
   *
   * Specializations of this template with D =1, 2, and 3 are derived 
   * from specializations of the base class template Rp::ScftThermo, and 
   * inherit their public interface and almost all of their source code
   * from this base class. See the documentation for this base class 
   * template for details. 
   *
   * \ingroup Rpc_Scft_Module
   */
   template <int D>
   class ScftThermo : public Rp::ScftThermo<D, System<D> >
   {
   public:

      /**
      * Constructor.
      *
      * \param system  parent System
      */
      ScftThermo(System<D> const & system);

   };

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::ScftThermo<1, Rpc::System<1> >;
      extern template class Rp::ScftThermo<2, Rpc::System<2> >;
      extern template class Rp::ScftThermo<3, Rpc::System<3> >;
   }
   namespace Rpc {
      extern template class ScftThermo<1>;
      extern template class ScftThermo<2>;
      extern template class ScftThermo<3>;
   }
} 
#endif
