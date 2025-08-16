#ifndef RPC_SCFT_THERMO_H
#define RPC_SCFT_THERMO_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/system/ScftThermoTmpl.h>    // base class template
#include <rpc/system/System.h>       // template parameter

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * Computes SCFT free energies.
   *
   * \ingroup Rpc_Scft_Module
   */
   template <int D>
   class ScftThermo : public ScftThermoTmpl<D, System<D> >
   {
   public:

      /// Alias for base class.
      using Base = ScftThermoTmpl<D, System<D> >;

      /**
      * Constructor.
      *
      * \param system  parent System
      */
      ScftThermo(System<D> const & system);

   protected:

      /// Alias for r-grid field type.
      using FieldT = typename Base::FieldT;
 
      /**
      * Inner product of fields (sum of elements on a grid).
      * 
      * \param A 1st field
      * \param B 2nd field
      */
      double innerProduct(FieldT const & A,
                          FieldT const & B) const override;

   };

   // Suppress implicit instantiation
   extern template class ScftThermo<1>;
   extern template class ScftThermo<2>;
   extern template class ScftThermo<3>;

} // namespace Rpc

namespace Prdc {
   // Suppress implicit instantiation of base class
   extern template class ScftThermoTmpl<1, Rpc::System<1> >;
   extern template class ScftThermoTmpl<2, Rpc::System<2> >;
   extern template class ScftThermoTmpl<3, Rpc::System<3> >;
} // namespace Rpc

} // namespace Pscf
#endif
