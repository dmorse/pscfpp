#ifndef RPC_EXPLICIT_BD_STEP_H
#define RPC_EXPLICIT_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/ExplicitBdStep.h>  // base class template
#include <rpc/system/Types.h>                // base template argument 
#include <rpc/fts/brownian/BdStep.h>         // indirect base class
#include <prdc/field/cpu/RField.h>                 // base class member

namespace Pscf {
namespace Rpc {

   // Forward declaration
   template <int D> class BdSimulator;

   /**
   * Explicit Euler-Maruyama Brownian dynamics step.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::ExplicBdStep, and
   * inherit their public interface and almost all of their source code
   * from this base class.  See the documentation of this base class 
   * template for details. 
   *
   * \see \ref rp_ExplicitBdStep_page "Manual Page"
   * \ingroup Rpc_Fts_Brownian_Module
   */
   template <int D>
   class ExplicitBdStep : public Rp::ExplicitBdStep<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator object
      */
      ExplicitBdStep(Rp::BdSimulator<D, Rpc::Types<D> >& simulator);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::ExplicitBdStep<1, Rpc::Types<1> >;
      extern template class Rp::ExplicitBdStep<2, Rpc::Types<2> >;
      extern template class Rp::ExplicitBdStep<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class ExplicitBdStep<1>;
      extern template class ExplicitBdStep<2>;
      extern template class ExplicitBdStep<3>;
   }
}
#endif
