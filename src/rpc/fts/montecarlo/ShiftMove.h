#ifndef RPC_SHIFT_MOVE_H
#define RPC_SHIFT_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/ShiftMove.h>     // base class template
#include <rpc/system/Types.h>                // base class argument
#include <rpc/fts/montecarlo/McMove.h>       // indirect base class
#include <prdc/cpu/RField.h>                 // base class member

namespace Pscf {
namespace Rpc {

   // Forward declarations
   template <int D> class McSimulator;
   template <int D> class System;

   using namespace Util;
   using namespace Prdc;

   /**
   * ShiftMove shifts field.
   *
   * \see \ref rp_ShiftMove_page "Manual Page".
   * \ingroup Rpc_Fts_MonteCarlo_Module
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

   protected:

      using McMove<D>::system;

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

      using RpShiftMove = Rp::ShiftMove<D, Types<D> >;

   };


}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::ShiftMove<1, Rpc::Types<1> >;
      extern template class Rp::ShiftMove<2, Rpc::Types<2> >;
      extern template class Rp::ShiftMove<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class ShiftMove<1>;
      extern template class ShiftMove<2>;
      extern template class ShiftMove<3>;
   }
}
#endif
