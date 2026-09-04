#ifndef RPC_SHIFT_MOVE_C_H
#define RPC_SHIFT_MOVE_C_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/ShiftMoveBase.h> // base class template
#include <pscf/backend/cpp/CPT.h>               // base class argument

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
   class ShiftMove<D,CPT>
     : public ShiftMoveBase<D,CPT>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator
      */
      ShiftMove(McSimulator<D,CPT>& simulator);

   protected:

      using McMove<D,CPT>::system;

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

      using ShiftMoveBaseT = ShiftMoveBase<D,CPT>;

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE_CPP(ShiftMove)
   //extern template class ShiftMove<1,CPT>;
   //extern template class ShiftMove<2,CPT>;
   //extern template class ShiftMove<3,CPT>;

} // namespace Rp
} // namespace Pscf
#endif
