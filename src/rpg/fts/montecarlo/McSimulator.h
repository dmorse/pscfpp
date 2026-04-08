#ifndef RPG_MC_SIMULATOR_H
#define RPG_MC_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/McSimulator.h>      // base class template
#include <rpg/system/Types.h>                   // template argument
#include <rpg/fts/simulator/Simulator.h>        // indirect base class
#include <rpg/fts/montecarlo/McMoveManager.h>   // member of base class
#include <rpg/fts/analyzer/AnalyzerManager.h>   // member of base class

namespace Pscf {
namespace Rpg {

   /**
   * Monte Carlo simulator for PS-FTS.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of base class template Rp::McSimulator, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see Rp::McSimulator
   * \see \ref rp_McSimulator_page (manual page)
   * \ingroup Rpg_Fts_MonteCarlo_Module
   */
   template <int D>
   class McSimulator : public Rp::McSimulator< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent system object
      */
      McSimulator(System<D>& system);

      McSimulator() = delete;
      McSimulator(McSimulator<D> const & ) = delete;

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class McSimulator<1, Rpg::Types<1> >;
      extern template class McSimulator<2, Rpg::Types<2> >;
      extern template class McSimulator<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class McSimulator<1>;
      extern template class McSimulator<2>;
      extern template class McSimulator<3>;
   }
}
#endif
