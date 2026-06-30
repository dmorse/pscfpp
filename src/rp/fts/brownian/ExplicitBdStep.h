#ifndef RP_EXPLICIT_BD_STEP_H
#define RP_EXPLICIT_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/BdStep.h>    // base class
#include <util/containers/DArray.h>    // member

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Explicit Euler-Maruyama Brownian dynamics step.
   *
   * Template parameters:
   *    - D : dimension
   *    - T : Types class, Rpc::Types<D> or Rpg::Types<D>
   *
   * \see \ref rp_ExplicitBdStep_page "Manual Page"
   * \ingroup Rp_Fts_Brownian_Module
   */
   template <int D, class T>
   class ExplicitBdStep : public BdStep<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator object
      */
      ExplicitBdStep(BdSimulator<D,T>& simulator);

      ~ExplicitBdStep() = default;

      // Not copyable or assignable
      ExplicitBdStep(ExplicitBdStep<D,T> const & ) = delete;
      ExplicitBdStep<D,T>& operator= (ExplicitBdStep<D,T> const &) = delete;

      /**
      * Read body of parameter file block.
      *
      * \param in  input parameter stream
      */
      void readParameters(std::istream &in) override;

      /**
      * Setup before simulation.
      */
      void setup() override;

      /**
      * Take a single Brownian dynamics step.
      *
      * \return true iff the compressor converged, false otherwise
      */
      bool step() override;

   protected:

      using BdStep<D,T>::system;
      using BdStep<D,T>::simulator;
      using BdStep<D,T>::vecRandom;

   private:

      // Local copy of w fields
      DArray< typename T::RField > w_;

      // Change in one component of wc
      typename T::RField dwc_;

      // Normal distributed random numbers
      typename T::RField gaussianField_;

      // Prefactor of -dc_ in deterministic drift term
      double mobility_;

   };

}
}
#endif
