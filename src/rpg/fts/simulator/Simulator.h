#ifndef RPG_SIMULATOR_H
#define RPG_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/simulator/Simulator.h>    // base class template
#include <rpg/system/Types.h>              // template argument
#include <rpg/fts/simulator/SimState.h>    // member
#include <prdc/cuda/RField.h>              // member (template arg)

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;

   /**
   * Field theoretic simulator (base class).
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::Simulator, and inherit
   * their entire public interface and almost all of their source code
   * from this base class.  
   *
   * \see Rp::Simulator
   * \ingroup Rpg_Fts_Simulator_Module
   */
   template <int D>
   class Simulator : public Pscf::Rp::Simulator<D, Types<D> >
   {
   public:

      /// Alias for direct base class.
      using Base = Pscf::Rp::Simulator<D, Types<D> >;

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      Simulator(Rp::System<D, Rpg::Types<D> >& system);

      /**
      * Destructor.
      */
      virtual ~Simulator() = default;

   protected:

      /**
      * Initialize seed for vector random number generator.
      */
      void initializeVecRandom();

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Simulator<1, Rpg::Types<1> >;
      extern template class Simulator<2, Rpg::Types<2> >;
      extern template class Simulator<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class Simulator<1>;
      extern template class Simulator<2>;
      extern template class Simulator<3>;
   }
}
#endif
