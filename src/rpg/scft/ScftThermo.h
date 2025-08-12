#ifndef RPG_SCFT_THERMO_H
#define RPG_SCFT_THERMO_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/system/ScftReal.h>   // base class template
#include <rpg/system/System.h>             // template parameter

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
   class ScftThermo : public ScftReal<D, System<D> >
   {
   public:

      /// Alias for base class
      using Base = ScftReal<D, System<D> >;

      /**
      * Constructor
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
      * \return sum of product A[i]*B[i] of values at mesh nodes.
      */
      double innerProduct(FieldT const & A,
                          FieldT const & B) const override;

   };

   // Suppress implicit instantiation
   extern template class ScftThermo<1>;
   extern template class ScftThermo<2>;
   extern template class ScftThermo<3>;

}

namespace Prdc {
   // Suppress implicit instantiation of base class
   extern template class ScftReal<1, Rpg::System<1> >;
   extern template class ScftReal<2, Rpg::System<2> >;
   extern template class ScftReal<3, Rpg::System<3> >;
}

} // namespace Pscf
#endif
