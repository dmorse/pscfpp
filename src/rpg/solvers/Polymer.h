#ifndef RPG_POLYMER_H
#define RPG_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Polymer.h>    // base class template
#include <rpg/system/Types.h>      // base class template parameter

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * Descriptor and solver for one polymer species.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of base class template Rp::Polymer, and
   * inherit their public interface and all of their source code from this
   * base class.
   *
   * \see Rp::Polymer
   * \ref user_param_polymer_sec "Manual Page"
   * \ingroup Rpg_Solver_Module
   */
   template <int D>
   class Polymer : public Rp::Polymer<D, Types<D> >
   {};

}
}

// Explicit instantiation declarations
namespace Pscf {
   extern template class PolymerTmpl< Rpg::Block<1>, Rpg::Propagator<1> >;
   extern template class PolymerTmpl< Rpg::Block<2>, Rpg::Propagator<2> >;
   extern template class PolymerTmpl< Rpg::Block<3>, Rpg::Propagator<3> >;
   namespace Rp {
      extern template class Polymer<1, Rpg::Types<1> >;
      extern template class Polymer<2, Rpg::Types<2> >;
      extern template class Polymer<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class Polymer<1>;
      extern template class Polymer<2>;
      extern template class Polymer<3>;
   }
}
#endif
