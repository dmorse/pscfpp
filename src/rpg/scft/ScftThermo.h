#ifndef RPG_SCFT_THERMO_H
#define RPG_SCFT_THERMO_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/system/ScftThermoTmpl.h>  // base class template
#include <rpg/system/System.h>           // base class template argument

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * Computes SCFT free energies.
   *
   * \ingroup Rpg_Field_Module
   */
   template <int D>
   class ScftThermo : public ScftThermoTmpl<D, System<D> >
   {
   public:

      /// Alias for base class
      using Base = ScftThermoTmpl<D, System<D> >;

      /**
      * Constructor
      *
      * \param system  parent System
      */
      ScftThermo(System<D> const & system);

   protected:

      /// Alias for r-grid field type.
      using RFieldT = typename Base::RFieldT;
 
      /**
      * Inner product of fields (sum of elements on a grid).
      * 
      * \param A 1st field
      * \param B 2nd field
      * \return sum of product A[i]*B[i] of values at mesh nodes.
      */
      double innerProduct(RFieldT const & A,
                          RFieldT const & B) const override;

   };

   // Explicit instantiation declarations
   extern template class ScftThermo<1>;
   extern template class ScftThermo<2>;
   extern template class ScftThermo<3>;

}

namespace Prdc {

   // Explicit instantiation declarations for base class
   extern template class ScftThermoTmpl<1, Rpg::System<1> >;
   extern template class ScftThermoTmpl<2, Rpg::System<2> >;
   extern template class ScftThermoTmpl<3, Rpg::System<3> >;

}

} // namespace Pscf
#endif
