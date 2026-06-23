#ifndef RPG_LINEAR_SWEEP_H
#define RPG_LINEAR_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/sweep/LinearSweep.h>      // direct base class template
#include <rpg/system/Types.h>               // base class argument
#include <rpg/scft/sweep/Sweep.h>           // indirect base class
#include <rpg/scft/sweep/SweepParameter.h>  // indirect base member

namespace Pscf {
namespace Rpg {

   // Forward declaration
   template <int D> class System;

   /**
   * Sweep in which parameters vary linearly with sweep variable s.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of base class template Rp::LinearSweep, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see Rp::LinearSweep
   * \see \ref scft_sweep_linear_sec "Manual page"
   * \ingroup Rpg_Scft_Sweep_Module
   */
   template <int D>
   class LinearSweep : public Rp::LinearSweep<D, Types<D> >
   {
   public:

      /**
      * Constructor.
      *
      * \param system  parent system object
      */
      LinearSweep(Rp::System<D, Rpg::Types<D> >& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class LinearSweep<1, Rpg::Types<1> >;
      extern template class LinearSweep<2, Rpg::Types<2> >;
      extern template class LinearSweep<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class LinearSweep<1>;
      extern template class LinearSweep<2>;
      extern template class LinearSweep<3>;
   }
}
#endif
