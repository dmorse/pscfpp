#ifndef RPC_FORCE_BIAS_MOVE_H
#define RPC_FORCE_BIAS_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/ForceBiasMoveBase.h>   // base class template
#include <pscf/cpu/CppTp.h>                      // base class argument
#include <rpc/fts/montecarlo/McMove.h>             // indirect base class
#include <prdc/field/cpu/RField.h>                 // base class member

namespace Pscf {
namespace Rp {

   template <int D, class T> class ForceBiasMove;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * ForceBiasMove attempts a Brownian dynamics move.
   *
   * This class implements a Monte Carlo move in which the unconstrained
   * attempted move is created by an explicit Euler-Maruyama Brownian
   * dynamics step. 
   *
   * Specializations of this template with D=1, 2, and 3 are derived
   * from specializations of base class template Rp::ForceBiasMove,
   * and inherit their public interface and almost all of their source
   * code from this base class.
   *
   * \see Rp::ForceBiasMove
   * \see \ref rp_ForceBiasMove_page "Manual Page"
   * \ingroup Rp_Fts_MonteCarlo_Module
   */
   template <int D>
   class ForceBiasMove<D, CppTp<D> > 
     : public Rp::ForceBiasMoveBase<D, CppTp<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator
      */
      ForceBiasMove(Rp::McSimulator<D, CppTp<D> >& simulator);

   private:

      /**
      * Compute force bias field.
      */
      void computeForceBias(RField<D, CppTp<D> >& result,
                            RField<D, CppTp<D> > const & di,
                            RField<D, CppTp<D> > const & df,
                            RField<D, CppTp<D> > const & dwc,
                            double mobility);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class ForceBiasMoveBase<1, CppTp<1> >;
      extern template class ForceBiasMoveBase<2, CppTp<2> >;
      extern template class ForceBiasMoveBase<3, CppTp<3> >;
      extern template class ForceBiasMove<1, CppTp<1> >;
      extern template class ForceBiasMove<2, CppTp<2> >;
      extern template class ForceBiasMove<3, CppTp<3> >;
   }
}
#endif
