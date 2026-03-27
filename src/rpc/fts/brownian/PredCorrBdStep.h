#ifndef RPC_PRED_CORR_BD_STEP_H
#define RPC_PRED_CORR_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/PredCorrBdStep.h> // base class template
#include <rpc/system/Types.h>               // base class template argument 
#include <prdc/cpu/RField.h>                // base class member
#include <rpc/fts/brownian/BdStep.h>        // indirect base class

namespace Pscf {
namespace Rpc {

   // Forward declaration
   template <int D> class BdSimulator;

   /**
   * Predictor-corrector Brownian dynamics time stepper.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::PredCorrBdStep, and
   * inherit their public interface and almost all of their source code
   * from this base class.  See the documentation of this base class 
   * template for details. 
   *
   * \see Rp::PredCorrBdStep
   * \see \ref rp_PredCorrBdStep_page "Manual Page"
   * \ingroup Rpc_Fts_Brownian_Module
   */
   template <int D>
   class PredCorrBdStep : public Rp::PredCorrBdStep<D, Types<D> >
   {
   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator object
      */
      PredCorrBdStep(BdSimulator<D>& simulator);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::PredCorrBdStep<1, Rpc::Types<1> >;
      extern template class Rp::PredCorrBdStep<2, Rpc::Types<2> >;
      extern template class Rp::PredCorrBdStep<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class PredCorrBdStep<1>;
      extern template class PredCorrBdStep<2>;
      extern template class PredCorrBdStep<3>;
   }
}
#endif
