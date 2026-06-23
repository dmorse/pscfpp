#ifndef RPG_BD_SIMULATOR_H
#define RPG_BD_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/BdSimulator.h>       // base class template
#include <rpg/system/Types.h>                  // template argument
#include <rpg/fts/simulator/Simulator.h>       // indirect base class
#include <rpg/fts/analyzer/AnalyzerManager.h>  // member of base class

namespace Pscf {
namespace Rpg {

   /**
   * Brownian dynamics simulator for PS-FTS.
   *
   * Specializations of this template with D=1, 2, and 3 are derived 
   * from specializations of the base class template Rp::BdSimulator, and
   * inherit their public interface and almost all of their source code
   * from this base class.  See the documentation of this base class 
   * for details. 
   *
   * \see Rp::BdSimulator
   * \see \ref rp_BdSimulator_page "Manual page"
   * \ingroup Rpg_Fts_Brownian_Module
   */
   template <int D>
   class BdSimulator : public Rp::BdSimulator<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent Rp::System<D, Rpg::Types<D> > object
      */
      BdSimulator(Rp::System<D, Rpg::Types<D> >& system);

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BdSimulator<1, Rpg::Types<1> >;
      extern template class BdSimulator<2, Rpg::Types<2> >;
      extern template class BdSimulator<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class BdSimulator<1>;
      extern template class BdSimulator<2>;
      extern template class BdSimulator<3>;
   }
}
#endif
