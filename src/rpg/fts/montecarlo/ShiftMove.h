#ifndef RPG_SHIFT_MOVE_H
#define RPG_SHIFT_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/ShiftMoveBase.h> // base class template
#include <pscf/backends/CUT.h>                // base class argument
#include <rpg/fts/montecarlo/McMove.h>       // indirect base class
#include <prdc/field/cuda/RField.h>          // base class member
#include <pscf/cuda/HostDArray.h>            // member
#include <pscf/cuda/cudaTypes.h>             // member

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * ShiftMove shifts field.
   *
   * \see \ref rp_ShiftMove_page "Manual Page".
   * \ingroup Rp_Fts_MonteCarlo_Module
   */
   template <int D>
   class ShiftMove<D,CUT>
    : public ShiftMoveBase<D,CUT>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator
      */
      ShiftMove(McSimulator<D,CUT>& simulator);

      /**
      * Setup before simulation.
      */
      void setup() override;

   protected:

      using McMove<D,CUT>::system;
    
      /**
      * Compute and store shifted w fields.
      *
      * On return, shifted values of fields obtained from system().w()
      * are stored in the w_ member array.
      *
      * \param shift  vector of integer shift values
      */ 
      void shiftFields(IntVec<D> const & shift) override;

   private:

      // Work space on CPU for unshifted field
      HostDArray<cudaReal> wOld_;

      // Work space on CPU for a shifted field
      HostDArray<cudaReal> wNew_;

      using ShiftMoveBaseT = ShiftMoveBase<D,CUT>;

   };


}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class ShiftMoveBase<1,CUT>;
      extern template class ShiftMoveBase<2,CUT>;
      extern template class ShiftMoveBase<3,CUT>;
      extern template class ShiftMove<1,CUT>;
      extern template class ShiftMove<2,CUT>;
      extern template class ShiftMove<3,CUT>;
   }
}
#endif
