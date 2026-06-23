#ifndef RPG_SWEEP_H
#define RPG_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/sweep/Sweep.h>             // base class template
#include <rpg/system/Types.h>                // base class argument
#include <rpg/scft/sweep/BasisFieldState.h>  // indirect base argument

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * Solve a sequence of SCFT problems along a line in parameter space.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of the base class template Rp::Sweep, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see Rp::Sweep
   * \see \ref scft_sweep_page "Manual page"
   * \ingroup Rpg_Scft_Sweep_Module
   */
   template <int D>
   class Sweep : public Rp::Sweep<D, Types<D> >
   {

   public:

      /**
      * Default constructor.
      */
      Sweep();

      /**
      * Constructor, creates assocation with parent system.
      *
      * \param system  parent system
      */
      Sweep(Rp::System<D, Rpg::Types<D> >& system);

      /**
      * Destructor.
      */
      virtual ~Sweep() = default;

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   extern template class SweepTmpl< Rpg::BasisFieldState<1> >;
   extern template class SweepTmpl< Rpg::BasisFieldState<2> >;
   extern template class SweepTmpl< Rpg::BasisFieldState<3> >;
   namespace Rp {
      extern template class Sweep<1, Rpg::Types<1> >;
      extern template class Sweep<2, Rpg::Types<2> >;
      extern template class Sweep<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class Sweep<1>;
      extern template class Sweep<2>;
      extern template class Sweep<3>;
   }
}
#endif
