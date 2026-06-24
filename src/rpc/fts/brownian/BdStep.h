#ifndef RPC_BD_STEP_H
#define RPC_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/BdStep.h>
#include <rpc/system/Types.h>

#if 0
namespace Pscf {
namespace Rpc {

   // Forward declaration
   template <int D> class BdSimulator;

   /**
   * BdStep is an abstract base class for Brownian dynamics steps.
   *
   * The virtual step() method must generate a single step.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::BdStep, and inherit
   * their entire public interface and almost all of their source code
   * from this base class.  
   *
   * \see \ref psfts_algo_brownian_page "Manual Page"
   * \ingroup Rpc_Fts_Brownian_Module
   */
   template <int D>
   class BdStep : public Rp::BdStep<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Rp::BdSimulator<D, Rpc::Types<D> > object
      */
      BdStep(Rp::BdSimulator<D, Rpc::Types<D> >& simulator);

      /**
      * Destructor.
      */
      virtual ~BdStep() = default;

   };

}
}
#endif

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BdStep< 1, Rpc::Types<1> >;
      extern template class BdStep< 2, Rpc::Types<2> >;
      extern template class BdStep< 3, Rpc::Types<3> >;
   }
}
#endif
