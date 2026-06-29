#ifndef RPC_LM_BD_STEP_H
#define RPC_LM_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/LMBdStep.h>     // base class template
#include <rpc/system/Types.h>             // base class template argument 
#include <prdc/field/cpu/RField.h>              // base class member
#include <rpc/fts/brownian/BdStep.h>      // indirect base class

namespace Pscf {
namespace Rpc {

   // Forward declaration
   template <int D> class BdSimulator;

   /**
   * Leimkuhler-Mathews Brownian dynamics time stepper.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::LMBdStep, and inherit
   * their entire public interface and almost all of their source code
   * from this base class.  See the documentation of this base class 
   * template for details. 
   *
   * \see Rp::LMBdStep
   * \see \ref rp_LMBdStep_page "Manual Page"
   * \ingroup Rpc_Fts_Brownian_Module
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
      LMBdStep(Rp::BdSimulator<D, Rpc::Types<D> >& simulator);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::LMBdStep<1, Rpc::Types<1> >;
      extern template class Rp::LMBdStep<2, Rpc::Types<2> >;
      extern template class Rp::LMBdStep<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class LMBdStep<1>;
      extern template class LMBdStep<2>;
      extern template class LMBdStep<3>;
   }
}
#endif
