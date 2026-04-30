#ifndef RPG_SCFT_THERMO_H
#define RPG_SCFT_THERMO_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/ScftThermo.h>         // base class template
#include <rpg/system/SystemConstRef.h>  // template argument
#include <rpg/system/Types.h>           // template argument

namespace Pscf {
namespace Rpg {

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
   * \ingroup Rpg_Scft_Module
   */
   template <int D>
   class ScftThermo : public Rp::ScftThermo<D, Types<D> >
   {
   public:

      /**
      * Constructor.
      *
      * \param system  parent System
      */
      ScftThermo(System<D> const & system);

      /**
      * Destructor.
      */
      virtual ~ScftThermo() = default;

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::ScftThermo<1, Rpg::Types<1> >;
      extern template class Rp::ScftThermo<2, Rpg::Types<2> >;
      extern template class Rp::ScftThermo<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class ScftThermo<1>;
      extern template class ScftThermo<2>;
      extern template class ScftThermo<3>;
   }
} 
#endif
