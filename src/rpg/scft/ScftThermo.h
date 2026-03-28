#ifndef RPG_SCFT_THERMO_H
#define RPG_SCFT_THERMO_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/ScftThermo.h>    // base class template
#include <rpg/system/System.h>     // base class template argument

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * Computes SCFT free energies.
   *
   * \ingroup Rpg_Field_Module
   */
   template <int D>
   class ScftThermo : public Rp::ScftThermo<D, System<D> >
   {
   public:

      /**
      * Constructor
      *
      * \param system  parent System
      */
      ScftThermo(System<D> const & system);

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class ScftThermo<1, Rpg::System<1> >;
      extern template class ScftThermo<2, Rpg::System<2> >;
      extern template class ScftThermo<3, Rpg::System<3> >;
   }
   namespace Rpg {
      extern template class ScftThermo<1>;
      extern template class ScftThermo<2>;
      extern template class ScftThermo<3>;
   }
}
#endif
