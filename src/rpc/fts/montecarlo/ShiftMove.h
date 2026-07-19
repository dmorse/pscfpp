#ifndef RPC_SHIFT_MOVE_H
#define RPC_SHIFT_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/ShiftMoveBase.h> // base class template
#include <pscf/cpu/Cpp.h>                // base class argument
#include <prdc/field/cpu/RField.h>           // base class member
#include <rpc/fts/montecarlo/McMove.h>       // indirect base class

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /**
   * ShiftMove shifts field.
   *
   * \see \ref rp_ShiftMove_page "Manual Page".
   * \ingroup Rp_Fts_MonteCarlo_Module
   */
   template <int D>
   class ShiftMove<D, Cpp<D> >
     : public ShiftMoveBase<D, Cpp<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator
      */
      ShiftMove(McSimulator<D, Cpp<D> >& simulator);

   protected:

      using McMove<D, Cpp<D> >::system;

   protected:
    
      /*
      * Compute and store shifted w fields.
      *
      * On return, shifted values of fields obtained from system().w()
      * are stored in the w_ member array.
      *
      * \param shift  vector of integer shift values
      */ 
      void shiftFields(IntVec<D> const & shift);

      using ShiftMoveBaseT = ShiftMoveBase<D, Cpp<D> >;

   };


}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class ShiftMoveBase<1, Cpp<1> >;
      extern template class ShiftMoveBase<2, Cpp<2> >;
      extern template class ShiftMoveBase<3, Cpp<3> >;
      extern template class ShiftMove<1, Cpp<1> >;
      extern template class ShiftMove<2, Cpp<2> >;
      extern template class ShiftMove<3, Cpp<3> >;
   }
}
#endif
