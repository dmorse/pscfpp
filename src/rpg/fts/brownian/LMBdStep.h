#ifndef RPG_LM_BD_STEP_H
#define RPG_LM_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/LMBdStep.h>     // base class template
#include <rpg/system/Types.h>             // base class template argument 
#include <rpg/fts/brownian/BdStep.h>      // indirect base class
#include <prdc/cuda/RField.h>              // base class member

namespace Pscf {
namespace Rpg {

   // Forward declaration
   template <int D> class BdSimulator;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * Leimkuhler-Mathews Brownian dynamics time stepper.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of base class template Rp::LMBdStep, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see Rp::LMBdStep
   * \see \ref rp_LMBdStep_page "Manual Page"
   * \ingroup Rpg_Fts_Brownian_Module
   */
   template <int D>
   class LMBdStep : public Rp::LMBdStep<D, Types<D> >
   {
   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator object
      */
      LMBdStep(BdSimulator<D>& simulator);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::LMBdStep<1, Rpg::Types<1> >;
      extern template class Rp::LMBdStep<2, Rpg::Types<2> >;
      extern template class Rp::LMBdStep<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class LMBdStep<1>;
      extern template class LMBdStep<2>;
      extern template class LMBdStep<3>;
   }
}
#endif
