#ifndef RPG_SHIFT_MOVE_H
#define RPG_SHIFT_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/ShiftMove.h>     // base class template
#include <rpg/system/Types.h>                // base class argument
#include <rpg/fts/montecarlo/McMove.h>       // indirect base class
#include <prdc/cuda/RField.h>                // base class member
#include <pscf/cuda/HostDArray.h>            // member
#include <pscf/cuda/cudaTypes.h>             // member

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;

   /**
   * ShiftMove shifts field.
   *
   * \see \ref rp_ShiftMove_page "Manual Page".
   * \ingroup Rpg_Fts_MonteCarlo_Module
   */
   template <int D>
   class ShiftMove : public Rp::ShiftMove<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator
      */
      ShiftMove(McSimulator<D>& simulator);

      /**
      * Setup before simulation.
      */
      void setup() override;

   protected:

      using McMove<D>::system;

    
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

      using RpShiftMove = Rp::ShiftMove<D, Types<D> >;

   };


}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::ShiftMove<1, Rpg::Types<1> >;
      extern template class Rp::ShiftMove<2, Rpg::Types<2> >;
      extern template class Rp::ShiftMove<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class ShiftMove<1>;
      extern template class ShiftMove<2>;
      extern template class ShiftMove<3>;
   }
}
#endif
