#ifndef RPG_SOLVENT_H
#define RPG_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Solvent.h>    // base class template
#include <rpg/system/Types.h>      // base class template parameter
#include <prdc/cuda/RField.h>      // member of base class

namespace Pscf {
namespace Rpg {

   /**
   * Solver and descriptor for a solvent species.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of the base class template Rp::Solvent,
   * and inherit their public interface and their entire implementation
   * from this base class.
   *
   * \see Rp::Solvent
   * \ref user_param_solvent_sec "Manual Page"
   * \ingroup Rpg_Solver_Module
   */
   template <int D>
   class Solvent : public Rp::Solvent<D, Rpg::Types<D> >
   {};

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Solvent<1, Rpg::Types<1> >;
      extern template class Solvent<2, Rpg::Types<2> >;
      extern template class Solvent<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class Solvent<1>;
      extern template class Solvent<2>;
      extern template class Solvent<3>;
   }
}
#endif
